"""
comparar_cineticas.py
Compara las cuatro cinéticas del proyecto en una muestra térmicamente delgada (0D, como un TGA):
la temperatura se impone (rampa + meseta) y solo cambia la química. Así se separa el efecto de la
cinética del efecto de la transferencia de calor en el reactor.

Casos de validación con datos experimentales de café (SCG):
  1. Inayat et al. (2022), Environ. Clim. Technol. 26(1) 720-729: lecho fijo de acero, café seco,
     ~0.9 °C/min hasta 300 °C, 30 min -> biochar 48.2 %, bio-oil ~22.5 %.
  2. Bok et al. (2012), Energy 47(1) 17-24: pirólisis rápida, máximo de bio-oil 54.85 % a 823 K.
  3. Programa por defecto del horno (10 °C/min, 500 °C, 30 min). Referencia genérica de pirólisis
     lenta (Inayat et al. 2022, citando [28]-[30]): ~30 % bio-oil, ~35 % biochar, ~35 % gas.

Ejecutar:  python comparar_cineticas.py
"""

import os

import numpy as np

from kinetics import PyrolysisKinetics
from kinetics_lumped import LumpedPyrolysisKinetics
from kinetics_rac import CineticaRAC

R_UNIVERSAL = 8.31446
T_AMB = 298.15

# Evaporación de humedad común a todos los modelos (misma ley que en el resto del proyecto)
A_EVAP, EA_EVAP = 5.0e5, 41.8e3


def k_evap(T):
    return A_EVAP * np.exp(-EA_EVAP / (R_UNIVERSAL * T)) if T > 350.0 else 0.0


def programa(beta_C_min, T_final_C, meseta_min):
    """Devuelve T(t) [K] y la duración total [s] de una rampa lineal + meseta."""
    beta = beta_C_min / 60.0
    T_f = T_final_C + 273.15
    t_rampa = (T_f - T_AMB) / beta
    t_total = t_rampa + meseta_min * 60.0
    return (lambda t: min(T_AMB + beta * t, T_f)), t_total, t_rampa


def dt_programa(beta_C_min):
    # Paso suficiente para seguir la rampa con < 0.5 K por paso
    return min(1.0, 0.5 / (beta_C_min / 60.0))


def simular_rac(regimen, T_de_t, t_total, dt, w_H2O, ceniza=0.0231):
    kin = CineticaRAC(regimen)
    seco = 1.0 - w_H2O
    m = kin.estado_inicial((1.0 - ceniza) * seco)    # base: 1 kg de café húmedo
    agua = w_H2O
    acum = {'gas': 0.0, 'agua_rxn': 0.0, 'bio_oil': 0.0}
    vapor = 0.0
    t = 0.0
    while t < t_total:
        h = min(dt, t_total - t)
        T = T_de_t(t)
        _, lib = kin.paso(m, T, h)
        for f in acum:
            acum[f] += float(lib[f])
        ev = agua * -np.expm1(-k_evap(T) * h)
        agua -= ev
        vapor += ev
        t += h
    sin_reaccionar = float(kin.masa_original(m))
    return {
        'biochar': float(kin.masa_solida(m)) - sin_reaccionar + ceniza * seco,
        'bio_oil': acum['bio_oil'],
        'agua': acum['agua_rxn'] + vapor,
        'gas': acum['gas'],
        'sin_reaccionar': sin_reaccionar,
        'humedad_restante': agua,
    }


def simular_interactivo_original(T_de_t, t_total, dt, w_H2O, beta_C_min, ceniza=0.0231):
    """Cinética 'calibrada' original de simulador_interactivo.py (sin fuente)."""
    A_c, E_c = 4.2e5, 88.0e3
    A_v, E_v = 2.8e10, 155.0e3
    A_g, E_g = 8.0e6, 105.0e3
    f_sec = np.clip(0.18 / np.sqrt(max(beta_C_min / 10.0, 0.2)), 0.04, 0.25)
    seco = 1.0 - w_H2O
    B, char, vol, gnc, agua, vapor = (1.0 - ceniza) * seco, ceniza * seco, 0.0, 0.0, w_H2O, 0.0
    t = 0.0
    while t < t_total:
        h = min(dt, t_total - t)
        T = T_de_t(t)
        kc, kv, kg = (A * np.exp(-E / (R_UNIVERSAL * T)) for A, E in ((A_c, E_c), (A_v, E_v), (A_g, E_g)))
        K = kc + kv + kg
        cons = B * -np.expm1(-K * h)
        B -= cons
        rv = cons * kv / K
        char += cons * kc / K + f_sec * rv
        vol += (1 - f_sec) * rv
        gnc += cons * kg / K
        ev = agua * -np.expm1(-k_evap(T) * h)
        agua -= ev
        vapor += ev
        t += h
    return {'biochar': char, 'bio_oil': vol, 'agua': vapor, 'gas': gnc, 'sin_reaccionar': B, 'humedad_restante': agua}


def simular_lumped(T_de_t, t_total, dt, w_H2O, ceniza=0.0231):
    """python/kinetics_lumped.py integrado con sub-pasos explícitos pequeños."""
    kin = LumpedPyrolysisKinetics()
    seco = 1.0 - w_H2O
    B, char, vol, gnc, agua, vapor = (1.0 - ceniza) * seco, ceniza * seco, 0.0, 0.0, w_H2O, 0.0
    t = 0.0
    while t < t_total:
        h = min(dt, t_total - t)
        T = T_de_t(t)
        r = kin.compute_rates(np.array(T), np.array(B), np.array(agua), np.array(0.0), np.array(char))
        k_tot = float(r['r_biomass_decay']) / max(B, 1e-30)
        cons = B * -np.expm1(-k_tot * h) if B > 0 else 0.0
        esc = cons / max(float(r['r_biomass_decay']) * h, 1e-30)   # escala las tasas a la masa consumida
        B -= cons
        char += float(r['d_char_dt']) * h * esc
        vol += float(r['d_volatiles_dt']) * h * esc
        gnc += float(r['d_gnc_dt']) * h * esc
        ev = agua * -np.expm1(-k_evap(T) * h)
        agua -= ev
        vapor += ev
        t += h
    return {'biochar': char, 'bio_oil': vol, 'agua': vapor, 'gas': gnc, 'sin_reaccionar': B, 'humedad_restante': agua}


def simular_multicomponente(T_de_t, t_total, dt, w_H2O):
    """python/kinetics.py: cada componente se reparte con rendimientos fijos + craqueo de alquitranes.
    Se usan los mismos parámetros de PyrolysisKinetics, integrados de forma exacta por componente."""
    from config import CoffeeBiomassConfig
    kin = PyrolysisKinetics()
    bio = CoffeeBiomassConfig(rho_bulk=1.0, moisture_fraction=w_H2O)
    comp = {
        'hce':  [bio.w_HCE,  kin.A_hce,  kin.Ea_hce,  kin.y_char_hce,  kin.y_tar_hce,  kin.y_gas_hce],
        'cell': [bio.w_CELL, kin.A_cell, kin.Ea_cell, kin.y_char_cell, kin.y_tar_cell, kin.y_gas_cell],
        'lig':  [bio.w_LIG,  kin.A_lig,  kin.Ea_lig,  kin.y_char_lig,  kin.y_tar_lig,  kin.y_gas_lig],
        'lip':  [bio.w_LIP,  kin.A_lip,  kin.Ea_lip,  kin.y_char_lip,  kin.y_tar_lip,  kin.y_gas_lip],
    }
    total0 = sum(c[0] for c in comp.values()) + bio.w_ASH + w_H2O   # la composición no suma 100 %
    char, tar, gas, agua, vapor = bio.w_ASH, 0.0, 0.0, w_H2O, 0.0
    t = 0.0
    while t < t_total:
        h = min(dt, t_total - t)
        T = T_de_t(t)
        RT = R_UNIVERSAL * T
        frac_crk = np.clip(kin.A_crk_gas * np.exp(-kin.Ea_crk_gas / RT) * 3.0, 0.0, 0.35)
        for c in comp.values():
            cons = c[0] * -np.expm1(-c[1] * np.exp(-c[2] / RT) * h)
            c[0] -= cons
            tar_p = c[4] * cons
            crk = frac_crk * tar_p
            char += c[3] * cons + 0.08 * crk
            tar += tar_p - crk - 0.08 * crk
            gas += c[5] * cons + crk
        ev = agua * -np.expm1(-k_evap(T) * h)
        agua -= ev
        vapor += ev
        t += h
    res = {'biochar': char, 'bio_oil': tar, 'agua': vapor, 'gas': gas,
           'sin_reaccionar': sum(c[0] for c in comp.values()), 'humedad_restante': agua}
    return {k: v / total0 for k, v in res.items()}


CASOS = [
    {
        'nombre': 'Inayat 2022: lecho fijo, 0.9 °C/min, 300 °C, 30 min (café seco)',
        'beta': 0.9, 'T_final': 300.0, 'meseta': 30.0, 'w_H2O': 0.0, 'regimen': 'lecho_fijo',
        'experimental': {'solido': 48.2, 'liquido': 22.5},
    },
    {
        'nombre': 'Bok 2012: pirólisis rápida a 550 °C (823 K)',
        'beta': 60000.0, 'T_final': 550.0, 'meseta': 1.0, 'w_H2O': 0.0, 'regimen': 'rapida',
        'experimental': {'liquido': 54.85},
    },
    {
        'nombre': 'Horno por defecto: 10 °C/min, 500 °C, 30 min (8 % humedad)',
        'beta': 10.0, 'T_final': 500.0, 'meseta': 30.0, 'w_H2O': 0.08, 'regimen': 'lecho_fijo',
        'experimental': {'solido': 35.0, 'liquido': 30.0, 'gas': 35.0},  # Referencia genérica, no de café
        'etiqueta_exp': 'Referencia\ngenérica (no café)',
    },
]


def correr_caso(caso):
    T_de_t, t_total, _ = programa(caso['beta'], caso['T_final'], caso['meseta'])
    dt = dt_programa(caso['beta'])
    w = caso['w_H2O']
    return {
        'RAC + TGL (nuevo)': simular_rac(caso['regimen'], T_de_t, t_total, dt, w),
        'Interactivo original': simular_interactivo_original(T_de_t, t_total, dt, w, caso['beta']),
        'Lumped': simular_lumped(T_de_t, t_total, min(dt, 0.05 if caso['beta'] > 100 else dt), w),
        'Multicomponente': simular_multicomponente(T_de_t, t_total, min(dt, 0.05 if caso['beta'] > 100 else dt), w),
    }


def tabla_markdown(caso, resultados):
    filas = ["| Modelo | Biochar | Sin reaccionar | **Sólido residual** | Bio-oil orgánico | Agua | **Líquido total** | **Gas** | Balance |",
             "|---|---|---|---|---|---|---|---|---|"]
    for modelo, r in resultados.items():
        liquido = r['bio_oil'] + r['agua']
        solido = r['biochar'] + r['sin_reaccionar']
        balance = sum(r.values())
        filas.append(f"| {modelo} | {100*r['biochar']:.1f} | {100*r['sin_reaccionar']:.1f} | **{100*solido:.1f}** | "
                     f"{100*r['bio_oil']:.1f} | {100*r['agua']:.1f} | **{100*liquido:.1f}** | **{100*r['gas']:.1f}** | {100*balance:.2f} |")
    exp = caso['experimental']
    filas.append(f"| *Experimental / referencia* | n/a | n/a | **{exp.get('solido', 'n/a')}** | n/a | n/a | "
                 f"**{exp.get('liquido', 'n/a')}** | **{exp.get('gas', 'n/a')}** | n/a |")
    return "\n".join(filas)


def graficar(todos, archivo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(todos), figsize=(6 * len(todos), 5.2), sharey=True)
    colores = {'solido': '#3E2723', 'liquido': '#BF360C', 'gas': '#1A237E'}
    for ax, (caso, resultados) in zip(axes, todos):
        modelos = list(resultados) + ['Experimental']
        etiquetas = [m.replace(' (nuevo)', '\n(nuevo)').replace(' original', '\noriginal') for m in modelos]
        etiquetas[-1] = caso.get('etiqueta_exp', 'Experimental')
        x = np.arange(len(modelos))
        ancho = 0.26
        for k, (clave, color) in enumerate(colores.items()):
            vals = []
            for m in modelos:
                if m == 'Experimental':
                    vals.append(caso['experimental'].get(clave, np.nan))
                else:
                    r = resultados[m]
                    if clave == 'liquido':
                        vals.append(100 * (r['bio_oil'] + r['agua']))
                    elif clave == 'solido':
                        vals.append(100 * (r['biochar'] + r['sin_reaccionar']))
                    else:
                        vals.append(100 * r[clave])
            ax.bar(x + (k - 1) * ancho, vals, ancho, color=color, edgecolor='black', label={'solido': 'Sólido residual', 'liquido': 'Líquido total', 'gas': 'Gas'}[clave])
        ax.set_xticks(x)
        ax.set_xticklabels(etiquetas, fontsize=8)
        ax.set_title(caso['nombre'], fontsize=9, fontweight='bold')
        ax.grid(axis='y', linestyle=':', alpha=0.6)
    axes[0].set_ylabel("Rendimiento (% en peso)")
    axes[0].legend()
    plt.tight_layout()
    plt.savefig(archivo, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    todos = []
    for caso in CASOS:
        res = correr_caso(caso)
        todos.append((caso, res))
        print(f"\n### {caso['nombre']}\n")
        print(tabla_markdown(caso, res))
    salida = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "figuras", "comparacion_cineticas.png")
    graficar(todos, salida)
    print(f"\nFigura: {salida}")

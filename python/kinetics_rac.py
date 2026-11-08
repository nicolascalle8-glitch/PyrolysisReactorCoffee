"""
kinetics_rac.py
Cinética detallada RAC (Ranzi - Anca-Couce, versión 2017) + triglicéridos (Debiagi et al. 2015)
para la pirólisis de residuos de café (SCG).

Fuentes de los parámetros (transcritos sin modificar):
  [RAC]  A. Anca-Couce (2017). "RAC detailed scheme for biomass pyrolysis", TU Graz. Reacciones R1-R25,
         tomadas de Anca-Couce, Sommersacher & Scharler (2017), J. Anal. Appl. Pyrolysis 127, 411-425
         y Anca-Couce & Scharler (2017), Fuel 206, 572-579. Base: Ranzi et al. (2008), Energy & Fuels 22, 4292.
  [TGL]  Debiagi et al. (2015), "Extractives Extend the Applicability of Multistep Kinetic Scheme of
         Biomass Pyrolysis", Energy & Fuels 29, 6544-6555, Tabla 5: TGL -> acroleína + 3 FFA.

Qué aporta frente a los modelos anteriores:
  - Reacciones competitivas con energías de activación distintas: la selectividad cambia con la
    temperatura y la velocidad de calentamiento (pirólisis lenta -> más char, rápida -> más líquido).
  - Carbonización secundaria heterogénea (parámetros x): alta en lechos fijos con partículas
    grandes y calentamiento lento (la canastilla), baja en pirólisis rápida.
  - Calor de reacción por reacción (kJ/g), positivo = endotérmico (convención del esquema RAC).

Integración: en cada paso cada especie se consume exactamente como B*(1 - e^(-K*dt)) y esa masa
se reparte entre los productos con rendimientos másicos que suman 1. La masa se conserva por
construcción para cualquier dt.
"""

import numpy as np

R_UNIVERSAL = 8.31446  # J/(mol·K)
MASA_ATOMICA = {'C': 12.011, 'H': 1.008, 'O': 15.999}

# Fórmulas (C, H, O) de las especies del esquema RAC y de Debiagi (2015)
FORMULAS = {
    # Sólidos
    'CELL': (6, 10, 5), 'CELLA': (6, 10, 5),
    'HCEA1': (5, 8, 4), 'HCEA2': (5, 8, 4),
    'LIG-C': (15, 14, 4), 'LIG-H': (22, 28, 9), 'LIG-O': (20, 22, 10),
    'LIG-CC': (15, 14, 4), 'LIG-OH': (19, 22, 8), 'LIG': (11, 12, 4),
    'G{CO2}': (1, 0, 2), 'G{CO}': (1, 0, 1), 'G{COH2}': (1, 2, 1), 'G2{COH2}': (1, 2, 1),
    'G3{COH2}': (1, 2, 1), 'G{H2}': (0, 2, 0), 'G{CH4}': (1, 4, 0), 'G{CH3OH}': (1, 4, 1),
    'G{C2H4}': (2, 4, 0), 'CHAR': (1, 0, 0),
    'TGL': (57, 100, 7),
    # Volátiles
    'HAA': (2, 4, 2), 'AA': (2, 4, 2), 'HCOOH': (1, 2, 2), 'GLYOX': (2, 2, 2), 'C3H6O': (3, 6, 1),
    'HMFU': (6, 6, 3), 'LVG': (6, 10, 5), 'XYL': (5, 8, 4), 'PCOUMARYL': (9, 10, 2),
    'PHENOL': (6, 6, 1), 'FE2MACR': (11, 12, 4), 'H2': (0, 2, 0), 'CO': (1, 0, 1), 'CO2': (1, 0, 2),
    'CH4': (1, 4, 0), 'CH2O': (1, 2, 1), 'CH3OH': (1, 4, 1), 'C2H4': (2, 4, 0), 'CH3CHO': (2, 4, 1),
    'ETOH': (2, 6, 1), 'H2O': (0, 2, 1), 'ACROLEINA': (3, 4, 1), 'FFA': (18, 32, 2),
}

SOLIDOS = ['CELL', 'CELLA', 'HCE', 'HCEA1', 'HCEA2', 'LIG-C', 'LIG-H', 'LIG-O', 'LIG-CC', 'LIG-OH', 'LIG',
           'G{CO2}', 'G{CO}', 'G{COH2}', 'G2{COH2}', 'G3{COH2}', 'G{H2}', 'G{CH4}', 'G{CH3OH}', 'G{C2H4}',
           'CHAR', 'TGL']
# Componentes de la biomasa original (lo que queda de ellos es "biomasa sin reaccionar")
ORIGINALES = ['CELL', 'HCE', 'LIG-C', 'LIG-H', 'LIG-O', 'TGL']

# Agrupación de los volátiles en las fracciones que se miden en el laboratorio
GASES = ['H2', 'CO', 'CO2', 'CH4', 'C2H4']                      # Gases no condensables
AGUA = ['H2O']                                                    # Agua de reacción
BIO_OIL = ['HAA', 'AA', 'HCOOH', 'GLYOX', 'C3H6O', 'HMFU', 'LVG', 'XYL', 'PCOUMARYL', 'PHENOL',
           'FE2MACR', 'CH2O', 'CH3OH', 'CH3CHO', 'ETOH', 'ACROLEINA', 'FFA']  # Orgánicos condensables
FRACCIONES = {'gas': GASES, 'agua_rxn': AGUA, 'bio_oil': BIO_OIL}

# Parámetros "x" de carbonización secundaria (tabla del esquema RAC 2017)
REGIMENES = {
    # Partículas y lecho grandes, calentamiento lento: canastilla en horno con rampa
    'lecho_fijo': {'CELL': 0.1, 'HCE': 0.2, 'LIG': 0.3, 'G': 0.4, 'azucares': False},
    # Pirólisis rápida: partículas finas, calentamiento de cientos de °C/s
    'rapida': {'CELL': 0.025, 'HCE': 0.05, 'LIG': 0.075, 'G': 0.1, 'azucares': True},
}

# Productos de carbonización secundaria (camino x) de cada familia
_X_CELL = {'CHAR': 5.5, 'H2O': 4, 'CO2': 0.5, 'H2': 1}
_X_HCE = {'CHAR': 4.5, 'H2O': 3, 'CO2': 0.5, 'H2': 1}
_X_LIGCC = {'CHAR': 15, 'H2O': 4, 'H2': 3}
_X_LIG = {'CHAR': 10.5, 'H2O': 3, 'CO2': 0.5, 'H2': 3}

_R2 = {'HAA': 0.45, 'GLYOX': 0.2, 'C3H6O': 0.3, 'HMFU': 0.25, 'H2': 0.05, 'CO': 0.31, 'CO2': 0.41,
       'CH2O': 0.4, 'CH3OH': 0.15, 'CH3CHO': 0.1, 'H2O': 0.83, 'HCOOH': 0.02, 'G{H2}': 0.05,
       'G{CH4}': 0.2, 'CHAR': 0.61}
_R6 = {'CO': 0.5, 'CO2': 0.5, 'CH4': 0.325, 'CH2O': 0.8, 'CH3OH': 0.1, 'C2H4': 0.25, 'ETOH': 0.125,
       'H2O': 0.025, 'HCOOH': 0.025, 'G{CO2}': 0.275, 'G{COH2}': 0.4, 'G{H2}': 0.125,
       'G{CH3OH}': 0.45, 'CHAR': 0.875}


def _reacciones(regimen, acetilos_hce):
    """
    Lista de reacciones con coeficientes MOLARES. Cada entrada:
      (nombre, reactivo, A, b, E [kJ/mol], productos principales, familia x, productos x, dh_principal, dh_x)
    k = A * T^b * exp(-E/RT). Calor de reacción [kJ/g de reactivo] = dh_principal*(1-x) + dh_x*x
    """
    azucares = regimen['azucares']
    r3_prod, r3_dh = ({'LVG': 1}, 0.53) if azucares else (_R2, 0.56)
    r8_prod, r8_dh = ({'XYL': 1}, 0.77) if azucares else (_R6, 0.25)
    return [
        ('R1',  'CELL',   4e13,    0, 188.37, {'CELLA': 1}, None, None, 0.0, 0.0),
        ('R2',  'CELLA',  2e6,     0, 80.0,   _R2, 'CELL', _X_CELL, 0.56, -1.47),
        ('R3',  'CELLA',  4,       1, 41.86,  r3_prod, 'CELL', _X_CELL, r3_dh, -1.47),
        ('R5',  'HCE',    1e10,    0, 129.77, {'AA': acetilos_hce, 'HCEA1': 0.58, 'HCEA2': 0.42}, None, None, 0.0, 0.0),
        ('R6',  'HCEA1',  1.2e9,   0, 125.58, _R6, 'HCE', _X_HCE, 0.25, -1.42),
        ('R7',  'HCEA1',  0.15,    1, 33.5,
         {'CO': 0.1, 'CO2': 0.8, 'CH2O': 0.3, 'H2O': 0.25, 'HCOOH': 0.05, 'G{CO2}': 0.15, 'G{CO}': 0.15,
          'G{COH2}': 1.2, 'G{H2}': 0.2, 'G{CH4}': 0.625, 'G{C2H4}': 0.375, 'CHAR': 0.875},
         'HCE', _X_HCE, -0.64, -1.42),
        ('R8',  'HCEA1',  3,       1, 46.05,  r8_prod, 'HCE', _X_HCE, r8_dh, -1.42),
        ('R9',  'HCEA2',  0.5e10,  0, 138.14,
         {'HAA': 0.2, 'CO': 0.175, 'CO2': 0.275, 'CH2O': 0.5, 'ETOH': 0.1, 'H2O': 0.2, 'HCOOH': 0.025,
          'G{CO2}': 0.4, 'G{COH2}': 0.925, 'G{CH4}': 0.25, 'G{CH3OH}': 0.3, 'G{C2H4}': 0.275, 'CHAR': 1.0},
         'HCE', _X_HCE, -0.14, -1.42),
        ('R10', 'LIG-C',  1.33e15, 0, 203.02,
         {'LIG-CC': 0.35, 'PCOUMARYL': 0.1, 'PHENOL': 0.08, 'CO': 0.32, 'CH2O': 0.3, 'H2O': 1.0,
          'G{COH2}': 0.7, 'G{CH4}': 0.495, 'G{C2H4}': 0.41, 'CHAR': 5.735}, None, None, -0.47, 0.0),
        ('R11', 'LIG-H',  0.67e13, 0, 156.97,
         {'LIG-OH': 1.0, 'HAA': 0.25, 'C3H6O': 0.5, 'G{C2H4}': 0.5}, None, None, 0.10, 0.0),
        ('R12', 'LIG-O',  0.33e9,  0, 106.74, {'LIG-OH': 1.0, 'CO2': 1.0}, None, None, -0.21, 0.0),
        ('R13', 'LIG-CC', 3e7,     0, 131.86,
         {'HAA': 0.35, 'PCOUMARYL': 0.3, 'PHENOL': 0.2, 'CO': 0.4, 'CH4': 0.65, 'C2H4': 0.6, 'H2O': 0.7,
          'G{CO}': 0.4, 'G{COH2}': 1.0, 'CHAR': 6.75}, 'LIG', _X_LIGCC, -0.09, -1.30),
        ('R14', 'LIG-OH', 1e8,     0, 125.58,
         {'LIG': 1.0, 'CO': 0.55, 'CO2': 0.05, 'CH4': 0.1, 'CH3OH': 0.6, 'H2O': 0.9, 'HCOOH': 0.05,
          'G{CO}': 0.6, 'G{COH2}': 0.85, 'G{H2}': 0.1, 'G{CH4}': 0.35, 'G{CH3OH}': 0.3, 'G{C2H4}': 0.2,
          'CHAR': 4.15}, None, None, -0.17, 0.0),
        ('R15', 'LIG',    4,       1, 50.2,   {'FE2MACR': 1.0}, 'LIG', _X_LIG, 0.95, -1.52),
        ('R16', 'LIG',    0.4e9,   0, 125.58,
         {'C3H6O': 0.2, 'CO': 1.0, 'CH4': 0.2, 'CH2O': 0.2, 'CH3OH': 0.4, 'CH3CHO': 0.2, 'H2O': 0.95,
          'HCOOH': 0.05, 'G{CO}': 0.45, 'G{COH2}': 0.5, 'G{CH4}': 0.4, 'G{C2H4}': 0.65, 'CHAR': 5.5},
         'LIG', _X_LIG, -0.35, -1.52),
        ('R17', 'LIG',    0.083,   1, 33.5,
         {'CO': 0.4, 'CH4': 0.2, 'CH2O': 0.4, 'H2O': 0.6, 'G{CO}': 0.2, 'G{COH2}': 2.0, 'G{CH4}': 0.4,
          'G{CH3OH}': 0.4, 'G{C2H4}': 0.5, 'CHAR': 6.0}, 'LIG', _X_LIG, -0.50, -1.52),
        ('R18', 'G{CO2}', 1e5,     0, 100.46, {'CO2': 1.0}, None, None, 0.0, 0.0),
        ('R19', 'G{CO}',  3e13,    0, 209.3,  {'CO': 1.0}, 'G', {'CHAR': 0.5, 'CO2': 0.5}, 0.0, -3.08),
        ('R20', 'G{COH2}', 1e6,    0, 100.46,
         {'G2{COH2}': 0.75, 'H2': 0.25, 'CO': 0.125, 'CO2': 0.0625, 'CHAR': 0.0625}, None, None, 0.31, 0.0),
        ('R21', 'G{H2}',  1e12,    0, 313.96, {'H2': 1.0}, None, None, 0.0, 0.0),
        ('R22', 'G{CH4}', 2e13,    0, 300.0,  {'CH4': 1.0}, None, None, 0.0, 0.0),
        ('R23', 'G{CH3OH}', 1.2e13, 0, 209.3, {'CH3OH': 1.0}, 'G', {'CHAR': 1.0, 'H2O': 1.0, 'H2': 1.0}, 0.0, -1.27),
        ('R24', 'G{C2H4}', 1e6,    0, 100.46, {'C2H4': 0.3, 'CH4': 0.7, 'CHAR': 0.7}, None, None, 0.46, 0.0),
        ('R25', 'G2{COH2}', 1.5e9, 0, 209.3,  {'G3{COH2}': 0.2, 'CO': 0.8, 'H2': 0.8}, None, None, 0.0, 0.0),
        # Debiagi et al. (2015), Tabla 5: 45 700 kcal/kmol = 191.2 kJ/mol. Calor de reacción no reportado -> 0.
        ('R26', 'TGL',    7e12,    0, 45.7 * 4.184, {'ACROLEINA': 1.0, 'FFA': 3.0}, None, None, 0.0, 0.0),
    ]


def masa_molar(formula):
    c, h, o = formula
    return c * MASA_ATOMICA['C'] + h * MASA_ATOMICA['H'] + o * MASA_ATOMICA['O']


# Composición del residuo de café (base seca) usada en el resto del proyecto (GIEM + FTIR, config.py)
COMPOSICION_SCG = {'celulosa': 0.12, 'hemicelulosa': 0.35, 'lignina': 0.28, 'lipidos': 0.15, 'cenizas': 0.0231}
# Reparto de la lignina en LIG-C / LIG-H / LIG-O. No hay análisis elemental del café: se usa el reparto
# de madera dura del esquema RAC (6.0 / 7.0 / 9.0 de 22.0). [Supuesto: medir C/H/O del café para ajustarlo]
REPARTO_LIGNINA_HW = {'LIG-C': 6.0 / 22.0, 'LIG-H': 7.0 / 22.0, 'LIG-O': 9.0 / 22.0}


class CineticaRAC:
    def __init__(self, regimen='lecho_fijo', acetilos_hce=0.1, reparto_lignina=None):
        """
        regimen: 'lecho_fijo' (carbonización secundaria alta: canastilla, rampa lenta) o 'rapida'.
        acetilos_hce: moles de ácido acético por unidad C5H8O4 de hemicelulosa (RAC: 0.1 madera blanda,
                      0.4 madera dura). Los galactomananos del café son poco acetilados -> 0.1 [Supuesto].
        """
        if regimen not in REGIMENES:
            raise ValueError(f"Régimen '{regimen}' no válido. Opciones: {list(REGIMENES)}")
        self.regimen = regimen
        x = REGIMENES[regimen]
        self.reparto_lignina = reparto_lignina or REPARTO_LIGNINA_HW

        formulas = dict(FORMULAS)
        # La hemicelulosa incluye sus grupos acetilo: C5H8O4 + a*C2H4O2
        formulas['HCE'] = (5 + 2 * acetilos_hce, 8 + 4 * acetilos_hce, 4 + 2 * acetilos_hce)
        self.formulas = formulas

        self.especies = SOLIDOS + [e for e in FORMULAS if e not in SOLIDOS]
        self.idx = {e: i for i, e in enumerate(self.especies)}
        self.reacciones = _reacciones(x, acetilos_hce)
        n_e, n_r = len(self.especies), len(self.reacciones)

        self.A = np.zeros(n_r)
        self.b = np.zeros(n_r)
        self.E = np.zeros(n_r)
        self.reactivo = np.zeros(n_r, dtype=int)
        self.dh = np.zeros(n_r)                 # J/kg de reactivo (positivo = endotérmico)
        self.Y = np.zeros((n_e, n_r))           # Rendimientos MÁSICOS: kg de especie / kg de reactivo
        self.error_atomico = {}                 # Desbalance atómico de cada reacción (para las pruebas)

        for j, (nombre, reac, A, b, E, prod, fam, prod_x, dh_p, dh_x) in enumerate(self.reacciones):
            xv = x[fam] if fam else 0.0
            molar = {e: (1.0 - xv) * nu for e, nu in prod.items()}
            if fam:
                for e, nu in prod_x.items():
                    molar[e] = molar.get(e, 0.0) + xv * nu

            atomos_reac = np.array(formulas[reac], dtype=float)
            atomos_prod = sum(nu * np.array(formulas[e], dtype=float) for e, nu in molar.items())
            self.error_atomico[nombre] = np.max(np.abs(atomos_prod - atomos_reac) / np.maximum(atomos_reac, 1.0))

            mw_reac = masa_molar(formulas[reac])
            for e, nu in molar.items():
                self.Y[self.idx[e], j] += nu * masa_molar(formulas[e]) / mw_reac
            # Cierre exacto: los coeficientes publicados vienen redondeados a 2-3 cifras
            self.Y[:, j] /= np.sum(self.Y[:, j])

            self.A[j], self.b[j], self.E[j] = A, b, E * 1e3
            self.reactivo[j] = self.idx[reac]
            self.dh[j] = (dh_p * (1.0 - xv) + dh_x * xv) * 1e6  # kJ/g -> J/kg

        # Matriz de incidencia reactivo-reacción para sumar las constantes por especie
        self.M_reac = np.zeros((n_e, n_r))
        self.M_reac[self.reactivo, np.arange(n_r)] = 1.0
        # Matrices de agrupación de volátiles en fracciones medibles
        self.P = {f: np.array([1.0 if e in lista else 0.0 for e in self.especies]) for f, lista in FRACCIONES.items()}
        self.es_solido = np.array([e in SOLIDOS for e in self.especies], dtype=float)
        self.es_original = np.array([e in ORIGINALES for e in self.especies], dtype=float)

    def estado_inicial(self, rho_organico, composicion=None, forma=()):
        """
        Masa inicial (kg/m^3) de cada especie a partir de la densidad de materia orgánica seca.
        La composición orgánica (celulosa + hemicelulosa + lignina + lípidos) se normaliza a 1:
        en el café queda ~7.7 % sin identificar (proteínas, extractivos), que se reparte en proporción
        a los componentes conocidos. [Supuesto: no hay cinética de proteínas en el esquema RAC]
        """
        c = composicion or COMPOSICION_SCG
        organico = c['celulosa'] + c['hemicelulosa'] + c['lignina'] + c['lipidos']
        fr = {
            'CELL': c['celulosa'] / organico,
            'HCE': c['hemicelulosa'] / organico,
            'TGL': c['lipidos'] / organico,
        }
        for lig, f in self.reparto_lignina.items():
            fr[lig] = c['lignina'] / organico * f
        m = np.zeros((len(self.especies),) + tuple(forma))
        for e, f in fr.items():
            m[self.idx[e]] = rho_organico * f
        return m

    def constantes(self, T):
        """k_j(T) para todas las reacciones. T puede ser escalar o arreglo. Devuelve (n_r, *T.shape)."""
        T = np.asarray(T, dtype=float)
        Tr = T[np.newaxis, ...]
        forma = (-1,) + (1,) * T.ndim
        return self.A.reshape(forma) * Tr ** self.b.reshape(forma) * np.exp(-self.E.reshape(forma) / (R_UNIVERSAL * Tr))

    def paso(self, m, T, dt):
        """
        Avanza las especies sólidas un paso dt a temperatura T (constante en el paso).
        m: arreglo (n_especies, ...) en kg/m^3; se modifica en su lugar.
        Devuelve:
          q_rxn: calor generado promedio en el paso (W/m^3, positivo = exotérmico)
          liberado: dict {fraccion: kg/m^3 liberados en el paso} para 'gas', 'agua_rxn', 'bio_oil'
        """
        k = self.constantes(T)                                   # (n_r, ...)
        K = np.tensordot(self.M_reac, k, axes=(1, 0))            # (n_e, ...) suma de k por reactivo
        consumido_esp = m * -np.expm1(-K * dt)                   # masa consumida de cada especie
        K_reac = K[self.reactivo]
        fraccion = np.divide(k, K_reac, out=np.zeros_like(k), where=K_reac > 0)
        consumido = consumido_esp[self.reactivo] * fraccion      # (n_r, ...) masa por reacción

        producido = np.tensordot(self.Y, consumido, axes=(1, 0))  # (n_e, ...)
        m -= consumido_esp
        m += producido * self.es_solido.reshape((-1,) + (1,) * (m.ndim - 1))
        # Los volátiles salen del lecho (arrastre con N2): no se acumulan en el sólido
        liberado = {f: np.tensordot(P, producido, axes=(0, 0)) for f, P in self.P.items()}

        q_rxn = -np.tensordot(self.dh, consumido, axes=(0, 0)) / dt
        return q_rxn, liberado

    def masa_solida(self, m):
        return np.tensordot(self.es_solido, m, axes=(0, 0))

    def masa_original(self, m):
        return np.tensordot(self.es_original, m, axes=(0, 0))

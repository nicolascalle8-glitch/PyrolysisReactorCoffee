"""
graficas.py — Gráficas de un Resultado (la app las reutiliza).

graficar(resultado, carpeta) -> rutas de: temperaturas.png (setpoint, horno, pared, superficie del lecho,
centro y termopar de control vs. t), masa.png (masa residual vs. t) y mapa_final.png (T r–z al final, con
la interfaz lecho–acero marcada).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def graficar(resultado, carpeta: str) -> list[str]:
    os.makedirs(carpeta, exist_ok=True)
    h = resultado.historia
    rutas = []

    fig, ax = plt.subplots(figsize=(8, 5))
    for clave, estilo, nombre in (("T_set", "k--", "Setpoint"), ("T_horno", "m-", "Horno"),
                                  ("T_pared", "C1-", "Pared externa del tubo"),
                                  ("T_lecho_superficie", "r-", "Superficie del lecho"),
                                  ("T_centro", "b-", "Centro del lecho"),
                                  ("T_termopar_control", "g:", "Termopar de control")):
        ax.plot(h["t"], h[clave], estilo, label=nombre)
    ax.set_xlabel("Tiempo (min)"); ax.set_ylabel("Temperatura (°C)"); ax.grid(alpha=0.4); ax.legend(fontsize=8)
    rutas.append(_guardar(fig, carpeta, "temperaturas.png"))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(h["t"], h["masa_residual"], "k-")
    ax.set_xlabel("Tiempo (min)"); ax.set_ylabel("Masa residual (% de la carga)"); ax.set_ylim(0, 105); ax.grid(alpha=0.4)
    rutas.append(_guardar(fig, carpeta, "masa.png"))

    t, T = resultado.mapas[-1]
    R, Z = np.meshgrid(resultado.r_mm, resultado.z_cm, indexing="ij")
    fig, ax = plt.subplots(figsize=(6, 6))
    cs = ax.contourf(R, Z, T, levels=30, cmap="inferno")
    ax.axvline(resultado.r_mm[resultado.n_lecho - 1], color="cyan", ls="--", label="Interfaz lecho–acero")
    fig.colorbar(cs, label="Temperatura (°C)")
    ax.set_xlabel("r (mm)"); ax.set_ylabel("z (cm)"); ax.set_title(f"Mapa de temperatura, t = {t:.1f} min")
    ax.legend(loc="lower left", fontsize=8)
    rutas.append(_guardar(fig, carpeta, "mapa_final.png"))
    return rutas


def _guardar(fig, carpeta, nombre):
    ruta = os.path.join(carpeta, nombre)
    fig.tight_layout()
    fig.savefig(ruta, dpi=120)
    plt.close(fig)
    return ruta

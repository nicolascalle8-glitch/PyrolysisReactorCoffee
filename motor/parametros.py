"""
parametros.py — Lee parametros/base.toml, aplica cambios y valida.

Recibe: la ruta del TOML y, opcional, un dict de cambios {"grupo.nombre": valor}.
Devuelve: un objeto Parametros (p["grupo.nombre"] da el valor).
Valida: que existan todas las claves del motor, la unidad esperada (dict UNIDADES), el rango [min, max],
las opciones de los parámetros de texto y que la composición del lecho sume 100 % (tolerancia 0.01).
Cualquier problema lanza ErrorParametro con el nombre del parámetro y el motivo: el motor no corre.
"""

import copy
import sys
from dataclasses import dataclass, field

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

ETIQUETAS = ("Dato", "Medido", "Norma", "Supuesto", "Calculado", "Por validar")
COMPOSICION = ("celulosa", "hemicelulosa", "lignina", "lipidos", "otros", "cenizas")
TOLERANCIA_COMPOSICION = 0.01  # puntos porcentuales

# Unidad esperada de cada parámetro que usa el motor (si el TOML trae otra, no corre)
UNIDADES = {
    "programa.T_inicial": "°C",
    "programa.rampa": "°C/min",
    "programa.T_final": "°C",
    "programa.meseta": "min",
    "control.modo_ssr": "-",
    "control.Kp": "1/K",
    "control.Ki": "1/(K·s)",
    "control.Kd": "s/K",
    "control.t_ciclo": "s",
    "control.t_muestreo": "s",
    "horno.modo": "-",
    "horno.tau_piedra": "s",
    "horno.P_max": "W",
    "horno.C_h": "J/K",
    "horno.UA_p": "W/K",
    "horno.h_c": "W/(m²·K)",
    "horno.emisividad": "-",
    "horno.T_amb": "°C",
    "geometria.R_lecho": "m",
    "geometria.e_acero": "m",
    "geometria.L": "m",
    "malla.Nr_lecho": "-",
    "malla.Nr_acero": "-",
    "malla.Nz": "-",
    "lecho.rho_aparente": "kg/m³",
    "lecho.humedad": "%",
    "lecho.cenizas": "%",
    "lecho.celulosa": "%",
    "lecho.hemicelulosa": "%",
    "lecho.lignina": "%",
    "lecho.lipidos": "%",
    "lecho.otros": "%",
    "lecho.R_contacto_extra": "m²·K/W",
    "acero.rho": "kg/m³",
    "acero.cp": "J/(kg·K)",
    "acero.k": "W/(m·K)",
    "cinetica.regimen": "-",
    "cinetica.acetilos_hce": "mol/mol",
    "termopares.control_r": "m",
    "termopares.control_z": "m",
    "termopares.tau": "s",
    "numerica.dt0": "s",
    "numerica.dt_max": "s",
    "numerica.dT_max": "K",
    "numerica.dX_max": "-",
    "numerica.t_registro": "s",}
ENTEROS = ("malla.Nr_lecho", "malla.Nr_acero", "malla.Nz")


class ErrorParametro(ValueError):
    """Parámetro mal definido: el mensaje dice cuál y por qué."""


@dataclass
class Parametro:
    valor: float | str
    unidad: str
    etiqueta: str
    fuente: str
    min: float | None = None
    max: float | None = None
    opciones: list = field(default_factory=list)


class Parametros:
    """Conjunto validado de parámetros. p["grupo.nombre"] -> valor; p.info(clave) -> Parametro."""

    def __init__(self, datos: dict[str, Parametro]):
        self._datos = datos
        _validar(self._datos)

    def __getitem__(self, clave: str):
        if clave not in self._datos:
            raise ErrorParametro(f"{clave}: el parámetro no existe")
        return self._datos[clave].valor

    def __contains__(self, clave):
        return clave in self._datos

    def claves(self) -> list[str]:
        return list(self._datos)

    def info(self, clave: str) -> Parametro:
        if clave not in self._datos:
            raise ErrorParametro(f"{clave}: el parámetro no existe")
        return self._datos[clave]

    def supuestos(self) -> list[str]:
        """Claves con etiqueta Supuesto o Por validar (la app las muestra)."""
        return [c for c, d in self._datos.items() if d.etiqueta in ("Supuesto", "Por validar")]

    def con_cambios(self, cambios: dict | None) -> "Parametros":
        """Copia con los valores cambiados; lanza ErrorParametro si una clave no existe o queda inválida."""
        datos = copy.deepcopy(self._datos)
        for clave, valor in (cambios or {}).items():
            if clave not in datos:
                raise ErrorParametro(f"{clave}: el parámetro no existe (¿error de tipeo?)")
            datos[clave].valor = valor
        return Parametros(datos)

    def __repr__(self):
        return f"Parametros({len(self._datos)} parámetros)"


def _validar(datos: dict[str, Parametro]) -> None:
    faltan = [c for c in UNIDADES if c not in datos]
    if faltan:
        raise ErrorParametro(f"{faltan[0]}: falta en el archivo de parámetros")
    for clave, d in datos.items():
        if clave in UNIDADES and d.unidad != UNIDADES[clave]:
            raise ErrorParametro(f"{clave}: unidad '{d.unidad}' no es la esperada '{UNIDADES[clave]}'")
        if d.etiqueta not in ETIQUETAS:
            raise ErrorParametro(f"{clave}: etiqueta '{d.etiqueta}' no válida; opciones {ETIQUETAS}")
        if d.opciones:
            if d.valor not in d.opciones:
                raise ErrorParametro(f"{clave}: '{d.valor}' no es una opción válida {d.opciones}")
            continue
        if isinstance(d.valor, bool) or not isinstance(d.valor, (int, float)):
            raise ErrorParametro(f"{clave}: el valor debe ser un número, no {d.valor!r}")
        if clave in ENTEROS:
            if float(d.valor) != int(d.valor):
                raise ErrorParametro(f"{clave}: debe ser entero, no {d.valor}")
            d.valor = int(d.valor)
        else:
            d.valor = float(d.valor)
        if d.valor != d.valor:  # NaN
            raise ErrorParametro(f"{clave}: el valor es NaN")
        if d.min is not None and d.valor < d.min:
            raise ErrorParametro(f"{clave}: {d.valor} {d.unidad} está por debajo del mínimo {d.min}")
        if d.max is not None and d.valor > d.max:
            raise ErrorParametro(f"{clave}: {d.valor} {d.unidad} está por encima del máximo {d.max}")
    suma = sum(datos[f"lecho.{c}"].valor for c in COMPOSICION)
    if abs(suma - 100.0) > TOLERANCIA_COMPOSICION:
        raise ErrorParametro(
            f"lecho: la composición (base seca) suma {suma:.3f} %, no 100 %. "
            "Ajusta lecho.otros = 100 − (celulosa + hemicelulosa + lignina + lipidos + cenizas)")


def cargar_parametros(ruta: str, cambios: dict | None = None) -> Parametros:
    """Lee el TOML de parámetros, aplica los cambios y valida todo."""
    with open(ruta, "rb") as f:
        crudo = tomllib.load(f)
    datos = {}
    for grupo, tabla in crudo.items():
        for nombre, d in tabla.items():
            clave = f"{grupo}.{nombre}"
            try:
                datos[clave] = Parametro(
                    valor=d["valor"], unidad=d["unidad"], etiqueta=d["etiqueta"], fuente=d["fuente"],
                    min=d.get("min"), max=d.get("max"), opciones=list(d.get("opciones", [])))
            except (KeyError, TypeError) as e:
                raise ErrorParametro(f"{clave}: falta el campo {e} (valor, unidad, etiqueta, fuente)") from None
    p = Parametros(datos)
    return p.con_cambios(cambios) if cambios else p

"""Motor del reactor vertical de pirólisis de café. Interfaz pública: ver spec §7."""
from .parametros import ErrorParametro, Parametros, cargar_parametros
from .resultados import Resultado
from .simulador import ErrorSimulacion, simular

__all__ = ["ErrorParametro", "ErrorSimulacion", "Parametros", "Resultado", "cargar_parametros", "simular"]

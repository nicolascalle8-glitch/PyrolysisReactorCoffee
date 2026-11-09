"""Motor del reactor vertical de pirólisis de café. Interfaz pública: ver spec §7."""
from .parametros import ErrorParametro, Parametros, cargar_parametros

__all__ = ["ErrorParametro", "Parametros", "cargar_parametros"]

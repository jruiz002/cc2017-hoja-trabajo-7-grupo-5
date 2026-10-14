"""Paquete compartido de la Hoja de Trabajo 7 (grupo 5): teoría de colas, simulador DES corregido y modelo SIR."""
from . import teoria, simulador, sir  # noqa: F401
from .teoria import mmc, mmcc, erlang_b, erlang_c, little, parametros_grupo  # noqa: F401
from .simulador import simular, replicar, resumir, ic_t, semillas_independientes  # noqa: F401

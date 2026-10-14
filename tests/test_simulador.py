"""Pruebas del simulador corregido. Los casos deterministas se calculan a mano (ver comentarios).

Ejecutar desde la raíz del repo:  .venv/bin/python -m pytest tests -q
"""
import math

import numpy as np
import pytest

from hospital_uci import simulador as S
from hospital_uci import teoria as T


def det(valor):
    """Gancho determinista: siempre devuelve ``valor`` (ignora el generador)."""
    return lambda rng: valor


# ---------------------------------------------------------------- casos deterministas (respuesta a mano)
def test_L_ponderado_por_tiempo_D_D_1():
    # Llegadas en t=2,4,...,100 (50 pacientes); cada uno ocupa la cama 0.5 días. Horizonte 101.
    # Cama ocupada 50*0.5 = 25 de 101 días => L = utilización = 25/101 (promedio en el tiempo).
    # Un promedio "por evento" daría (1+0)/2 = 0.5, muy distinto: esta prueba detecta el punto (i).
    r = S.simular(1, 1, 1, 101.0, 0.0, semilla=0, interarribo=det(2.0), servicio=det(0.5))
    assert r["L"] == pytest.approx(25 / 101)
    assert r["utilizacion"] == pytest.approx(25 / 101)
    assert r["Lq"] == 0 and r["Wq"] == 0 and r["p_esperar"] == 0
    assert r["n_llegadas"] == 50 and r["W"] == pytest.approx(0.5)


def test_calentamiento_descarta_cohorte_inicial():
    # Mismo caso con calentamiento 11: llegadas en 12,14,...,100 (45). Cama ocupada 22.5 de T=90 => 0.25.
    r = S.simular(1, 1, 1, 101.0, 11.0, semilla=0, interarribo=det(2.0), servicio=det(0.5))
    assert r["n_llegadas"] == 45
    assert r["L"] == pytest.approx(22.5 / 90)


def test_espera_desde_la_llegada_y_cola_D_D_1():
    # Llegadas en t=1,2,...,10; servicio 2 días; 1 cama; horizonte 10.5.
    # Paciente k empieza servicio en 2k-1 => espera k-1 (medida desde la llegada). Wq = (0+1+...+9)/10 = 4.5; W = 6.5.
    # Si la espera se midiera desde el inicio del servicio daría 0 (punto vii).
    # P(esperar) = 9/10 (todos menos el primero).
    r = S.simular(1, 1, 1, 10.5, 0.0, semilla=0, interarribo=det(1.0), servicio=det(2.0))
    assert r["Wq"] == pytest.approx(4.5)
    assert r["W"] == pytest.approx(6.5)
    assert r["p_esperar"] == pytest.approx(0.9)
    # Área de cola a mano: pacientes 2..5 esperan 1+2+3+4 = 10; pacientes 6..10 esperan (10.5-k) hasta el horizonte: 12.5.
    assert r["Lq"] == pytest.approx(22.5 / 10.5)
    assert r["L"] == pytest.approx(32 / 10.5)       # cola 22.5 + cama ocupada 9.5 desde t=1
    assert r["n_dentro_al_horizonte"] == 6          # pacientes 5..10


def test_pacientes_al_final_se_drenan_sin_censura():
    # Mismo caso: solo 5 de los 10 pacientes han empezado servicio al llegar el horizonte.
    # Descartar los no atendidos daría Wq = (0+1+2+3+4)/5 = 2 (sesgo a la baja); el tratamiento consistente da 4.5.
    r = S.simular(1, 1, 1, 10.5, 0.0, semilla=0, interarribo=det(1.0), servicio=det(2.0))
    assert r["Wq"] != pytest.approx(2.0)
    assert r["n_admitidos"] == 10


def test_bloqueo_D_D_1_1():
    # 1 cama, sin espera, llegadas en 1,2,...,10, servicio 1.5, horizonte 10.2.
    # Se admiten los pacientes de t=1,3,5,7,9 y se bloquean 2,4,6,8,10 => P(bloqueo)=0.5; tasa efectiva = 5/10.2.
    # Ocupación: 4*1.5 + (10.2-9) = 7.2 => camas ocupadas = 7.2/10.2; W = 1.5.
    r = S.simular(1, 1, 1, 10.2, 0.0, semilla=0, con_espera=False, interarribo=det(1.0), servicio=det(1.5))
    assert r["p_bloqueo"] == pytest.approx(0.5)
    assert r["lam_efectiva"] == pytest.approx(5 / 10.2)
    assert r["camas_ocupadas"] == pytest.approx(7.2 / 10.2)
    assert r["W"] == pytest.approx(1.5)
    assert r["Lq"] == 0 and r["p_esperar"] == 0


def test_desempate_salida_antes_que_llegada():
    # Llegadas cada 2 días, servicio exactamente 2: la cama se libera en el instante en que llega el siguiente.
    # Con 1 cama y sin espera no debe haber bloqueos si la salida se procesa primero.
    r = S.simular(1, 1, 1, 21.0, 0.0, semilla=0, con_espera=False, interarribo=det(2.0), servicio=det(2.0))
    assert r["p_bloqueo"] == 0.0


# ---------------------------------------------------------------- reproducibilidad y semillas
def test_misma_semilla_mismo_resultado():
    a = S.simular(3.8, 0.25, 18, 300, 30, semilla=7)
    b = S.simular(3.8, 0.25, 18, 300, 30, semilla=7)
    assert a == b


def test_replicas_semillas_distintas_y_resultados_distintos():
    reps = S.replicar(30, 123, lam=3.8, mu=0.25, c=18, horizonte=200, calentamiento=20)
    assert len({r["semilla_estado"] for r in reps}) == 30
    assert len({round(r["L"], 12) for r in reps}) == 30


def test_replicas_no_correlacionadas_con_semillas_consecutivas():
    # Con semillas hijas de SeedSequence, el L de réplicas vecinas no debe estar correlacionado (|r| pequeño).
    reps = S.replicar(200, 99, lam=3.8, mu=0.25, c=18, horizonte=100, calentamiento=10)
    x = np.array([r["L"] for r in reps])
    r1 = np.corrcoef(x[:-1], x[1:])[0, 1]
    assert abs(r1) < 0.3  # 3.5 desv. típicas ≈ 3.5/sqrt(200) = 0.25


# ---------------------------------------------------------------- IC: unidad de muestreo = réplica
def test_ic_usa_n_replicas_no_n_pacientes():
    reps = S.replicar(30, 5, lam=3.8, mu=0.25, c=18, horizonte=300, calentamiento=30)
    res = S.resumir(reps)
    x = np.array([r["Wq"] for r in reps])
    h = 2.045229642 * x.std(ddof=1) / math.sqrt(30)  # t_{0.975, 29}
    assert res["Wq"]["n"] == 30
    assert res["Wq"]["semiancho"] == pytest.approx(h, rel=1e-6)


# ---------------------------------------------------------------- respuestas teóricas conocidas (estadísticas)
def _dentro_3se(reps, clave, teo):
    x = np.array([r[clave] for r in reps])
    se = x.std(ddof=1) / math.sqrt(len(x))
    return abs(x.mean() - teo) < 3 * se


def test_mm1_contra_formulas_cerradas():
    # M/M/1 con λ=0.5, μ=1: ρ=0.5, L=1, Lq=0.5, Wq=1, W=2, P(esperar)=0.5.
    reps = S.replicar(30, 2024, lam=0.5, mu=1.0, c=1, horizonte=6000, calentamiento=200)
    for clave, teo in {"L": 1.0, "Lq": 0.5, "Wq": 1.0, "W": 2.0, "p_esperar": 0.5, "utilizacion": 0.5}.items():
        assert _dentro_3se(reps, clave, teo), clave


def test_erlang_b_c1():
    # c=1, a=1 => B = a/(1+a) = 0.5;  L = a(1-B) = 0.5;  λ_ef = 0.5.
    reps = S.replicar(30, 2025, lam=1.0, mu=1.0, c=1, horizonte=6000, calentamiento=200, con_espera=False)
    for clave, teo in {"p_bloqueo": 0.5, "L": 0.5, "lam_efectiva": 0.5, "W": 1.0}.items():
        assert _dentro_3se(reps, clave, teo), clave


def test_validaciones_de_entrada():
    with pytest.raises(ValueError):
        S.simular(1, 1, 0, 10, 0)
    with pytest.raises(ValueError):
        S.simular(1, 1, 1, 10, 10)

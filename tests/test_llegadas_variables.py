# Pruebas con respuesta conocida del módulo llegadas_variables (adelgazamiento y cierres, Tasks 3.1c y 4.1b-c).
import math

import numpy as np
import pytest

from hospital_uci import llegadas_variables as LV
from hospital_uci import simulador as S
from hospital_uci import teoria as T


# Gancho determinista: siempre devuelve `valor` (ignora el generador).
def det(valor):
    return lambda rng: valor


# Caso a mano sin cierre: 2 camas, estancia 3, llegadas en 1, 2, 3.5, 4 y 4.5.
# t=1 y t=2 entran (salen en 4 y 5); t=3.5 se bloquea (2 ocupadas); t=4: sale el de t=1 antes (empate) y entra;
# t=4.5: ocupadas 2 -> bloqueado. Bloqueados = [F, F, T, F, T].
def test_caso_a_mano_sin_cierre():
    r = LV.simular_perdida(2, 20.0, semilla=0, servicio=det(3.0), llegadas=[1, 2, 3.5, 4, 4.5])
    assert r["bloqueado"].tolist() == [False, False, True, False, True]


# Caso a mano con cierre: 3 camas, estancia 10, llegadas en 0, 0.5 y 1; cierre de 2 camas en [2, 5).
# Nadie es expulsado: las 3 camas siguen ocupadas, así que t=3 (capacidad 1) y t=6 (capacidad 3, llenas) se bloquean;
# t=12 entra (ya salieron todos). En [2, 5) el hospital está "lleno" todo el tiempo y la ocupación media es 3.
def test_cierre_no_expulsa_y_bloquea():
    r = LV.simular_perdida(3, 30.0, semilla=0, servicio=det(10.0), llegadas=[0, 0.5, 1, 3, 6, 12],
                           cierres=[(2.0, 5.0, 2)])
    assert r["bloqueado"].tolist() == [False, False, False, True, True, False]
    lleno, ocup = LV.tiempo_lleno_por_ventana(r, [2.0, 5.0, 10.0, 11.0])
    assert lleno[0] == pytest.approx(1.0) and ocup[0] == pytest.approx(3.0)
    # en [5, 10) capacidad 3 y ocupadas 3 -> lleno; en [10, 11): salen en 10 y 10.5 -> ocupación media 1.5
    assert lleno[1] == pytest.approx(1.0)
    assert ocup[2] == pytest.approx((2 * 0.5 + 1 * 0.5) / 1.0)


# Durante el cierre la capacidad baja a c − k_off solo en [t0, t1).
def test_capacidad_en():
    cierres = [(10.0, 17.0, 6)]
    assert LV.capacidad_en(9.99, 23, cierres) == 23
    assert LV.capacidad_en(10.0, 23, cierres) == 17
    assert LV.capacidad_en(17.0, 23, cierres) == 23


# Con λ(t) escalonada (2 en [0, 50) y 6 en [50, 100)), el número esperado de llegadas es la integral de λ(t):
# 2·50 + 6·50 = 400 (100 y 300 por tramo). Se comprueba con 400 réplicas y un margen de 4 errores estándar.
def test_adelgazamiento_integral_escalonada():
    tasa = lambda t: 2.0 if t < 50 else 6.0  # noqa: E731
    hijas, _ = S.semillas_independientes(910001, 400)
    n1, n2 = [], []
    for h in hijas:
        t = LV.llegadas_adelgazamiento(np.random.default_rng(h), tasa, 6.0, 0.0, 100.0)
        n1.append(np.sum(t < 50)); n2.append(np.sum(t >= 50))
    for obs, esp in [(n1, 100.0), (n2, 300.0)]:
        ee = math.sqrt(esp / len(obs))          # Poisson: varianza = media
        assert abs(np.mean(obs) - esp) < 4 * ee
    assert np.var(n2, ddof=1) == pytest.approx(300.0, rel=0.25)   # dispersión de Poisson


# λ(t) > λmax debe detenerse con error (un λmax insuficiente sesgaría las llegadas a la baja sin avisar).
def test_adelgazamiento_rechaza_lam_max_insuficiente():
    with pytest.raises(ValueError):
        LV.llegadas_adelgazamiento(np.random.default_rng(1), LV.tasa_constante(5.0), 4.0, 0.0, 100.0)


# Con λ(t) constante = 3.8 generada por adelgazamiento (λmax = 5.2) y c = 18 se recupera Erlang B = 0.0911,
# tanto por la fracción de llegadas bloqueadas como por la fracción del tiempo con las camas llenas (PASTA).
def test_adelgazamiento_constante_reproduce_erlang_b():
    def resumen(r):
        m = r["t_llegada"] >= 100.0
        lleno, _ = LV.tiempo_lleno_por_ventana(r, [100.0, 3100.0])
        return {"p": r["bloqueado"][m].mean(), "lleno": lleno[0]}
    reps = LV.replicar_perdida(20, 910002, resumen, c=18, horizonte=3100.0,
                               tasa=LV.tasa_constante(3.8), lam_max=5.2)
    B = T.erlang_b(18, 15.2)
    for k in ("p", "lleno"):
        ic = S.ic_t([r[k] for r in reps])
        assert ic["ic_inf"] <= B <= ic["ic_sup"]


# La lognormal tiene media 1/μ = 4 y ln(estancia) tiene media m y desviación s (fórmulas del 2.2a).
def test_estancia_lognormal_momentos():
    g = LV.estancia_lognormal(0.25, 2.0)
    rng = np.random.default_rng(910003)
    x = np.array([g(rng) for _ in range(200_000)])
    s2 = math.log(5.0)
    assert np.mean(np.log(x)) == pytest.approx(math.log(4.0) - s2 / 2, abs=0.01)
    assert np.std(np.log(x)) == pytest.approx(math.sqrt(s2), abs=0.01)
    assert np.mean(x) == pytest.approx(4.0, rel=0.03)


# La tasa de la ola vale λ antes de que empiece la epidemia y λ + λepi dentro de la curva.
def test_tasa_ola():
    curva = {"t_dias": np.array([0.0, 1.0, 2.0]), "lambda_epi_pac_dia": np.array([0.0, 1.0, 0.5])}
    f = LV.tasa_ola(curva, 3.8, t_inicio=10.0)
    assert f(5.0) == 3.8
    assert f(10.5) == pytest.approx(4.3)
    assert f(11.5) == pytest.approx(4.55)


# Un "cierre" con k_off negativo agrega camas temporales (camas de contingencia del 4.1c) solo en esa ventana.
# 1 cama, estancia 10, llegadas en 0, 1 y 6; 1 cama extra en [5, 8): t=1 se bloquea y t=6 entra en la cama temporal.
def test_camas_temporales_con_k_off_negativo():
    r = LV.simular_perdida(1, 30.0, semilla=0, servicio=det(10.0), llegadas=[0, 1, 6], cierres=[(5.0, 8.0, -1)])
    assert LV.capacidad_en(6.0, 1, [(5.0, 8.0, -1)]) == 2
    assert r["bloqueado"].tolist() == [False, True, False]

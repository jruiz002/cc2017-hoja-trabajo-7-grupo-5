import math

import numpy as np
import pytest

from hospital_uci import teoria as T
from hospital_uci import sir


def test_parametros_grupo5():
    p = T.parametros_grupo(5)
    assert p["lam"] == pytest.approx(3.8) and p["c"] == 18 and p["a"] == pytest.approx(15.2)
    assert p["rho"] == pytest.approx(0.844444, abs=1e-6) and p["k_off"] == 6 and p["beta"] == pytest.approx(0.30)


def test_erlang_b_valores_conocidos():
    assert T.erlang_b(1, 1.0) == pytest.approx(0.5)
    assert T.erlang_b(2, 1.0) == pytest.approx(0.2)   # (1/2)/(1+1+1/2)
    assert T.erlang_b(0, 3.0) == 1.0


@pytest.mark.parametrize("c,a", [(18, 15.2), (5, 3.0), (30, 25.0), (1, 0.4)])
def test_recurrencia_igual_a_formula_directa(c, a):
    assert T.erlang_b(c, a) == pytest.approx(T.erlang_b_directa(c, a), rel=1e-12)


def test_mm1_formulas():
    m = T.mmc(0.5, 1.0, 1)
    assert m["p_esperar"] == pytest.approx(0.5) and m["Lq"] == pytest.approx(0.5)
    assert m["Wq"] == pytest.approx(1.0) and m["W"] == pytest.approx(2.0) and m["L"] == pytest.approx(1.0)


def test_little_ambas_variantes_grupo5():
    m = T.mmc(3.8, 0.25, 18)
    assert T.little(m["L"], 3.8, m["W"])["cumple"]
    assert T.little(m["Lq"], 3.8, m["Wq"])["cumple"]
    b = T.mmcc(3.8, 0.25, 18)
    assert T.little(b["L"], b["lam_efectiva"], b["W"])["cumple"]
    # usar λ en lugar de λ_ef: λW = L/(1-B), error relativo = B/(1-B)
    B = b["p_bloqueo"]
    assert T.little(b["L"], 3.8, b["W"])["error_relativo"] == pytest.approx(B / (1 - B), rel=1e-9)


def test_sir_conserva_poblacion_y_converge():
    c1 = sir.integrar_sir(dt=0.1)
    assert np.allclose(c1["S"] + c1["I"] + c1["R"], sir.N_DEFECTO, atol=1e-6)
    ref = sir.integrar_referencia(dt=0.1)
    assert np.max(np.abs(c1["lambda_epi"] - ref["lambda_epi"])) < 1e-6
    p_a, p_b = sir.pico(sir.integrar_sir(dt=0.1)), sir.pico(sir.integrar_sir(dt=0.01))
    assert abs(p_a["valor"] - p_b["valor"]) / p_b["valor"] < 1e-4


def test_sir_solucion_cerrada_tamano_final():
    # Relación exacta del SIR: ln(S(t)/S(0)) = -(β/γ)·R(t)/N, válida en todo t (se evalúa a los 400 días).
    c = sir.integrar_sir(dt=0.1, dias=400)
    Sf = c["S"][-1]
    lhs = math.log(Sf / c["S"][0])
    rhs = -(0.30 / 0.2) * (c["R"][-1] / sir.N_DEFECTO)
    assert lhs == pytest.approx(rhs, rel=1e-4)

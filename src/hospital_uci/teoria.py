"""Teoría de colas: M/M/c (con espera) y M/M/c/c (sin espera).

Unidades: tasas en 1/día, tiempos devueltos en DÍAS salvo las claves terminadas en ``_h`` (horas).

Notación (diapositivas de la sesión):
    lam  = λ   tasa de llegadas (pacientes/día)
    mu   = μ   tasa de servicio por cama (1/día); estancia media = 1/μ
    c          número de camas (servidores)
    a    = λ/μ carga ofrecida (Erlangs = camas ocupadas en promedio si nadie se perdiera)
    rho  = a/c utilización por cama (solo tiene sentido en M/M/c y exige a < c)
"""
from __future__ import annotations

from math import isclose


def parametros_grupo(g: int) -> dict:
    """Parámetros del grupo g según la hoja (λ en pacientes/día, c en camas, β en 1/día)."""
    import math

    mu = 0.25
    lam = 3.0 + 0.2 * ((g - 1) % 5)
    c = 17 + ((g - 1) % 3)
    a = lam / mu
    rho = a / c
    k_off = math.ceil(c - a) + 3  # menor entero >= c - a, más 3
    beta = 0.30 + 0.02 * ((g - 1) % 4)
    return {"g": g, "lam": lam, "mu": mu, "c": c, "a": a, "rho": rho, "k_off": k_off, "beta": beta}


def carga_ofrecida(lam: float, mu: float) -> float:
    """a = λ/μ  (adimensional, 'Erlangs')."""
    return lam / mu


def erlang_b(c: int, a: float) -> float:
    """Probabilidad de bloqueo de Erlang B por recurrencia (estable numéricamente).

    B(0, a) = 1;   B(k, a) = a·B(k-1, a) / (k + a·B(k-1, a)),  k = 1..c
    """
    if c < 0 or a < 0:
        raise ValueError("c y a deben ser >= 0")
    b = 1.0
    for k in range(1, int(c) + 1):
        b = a * b / (k + a * b)
    return b


def erlang_b_directa(c: int, a: float) -> float:
    """B(c, a) = (a^c / c!) / Σ_{k=0}^{c} a^k / k!  — solo para contrastar la recurrencia."""
    from math import factorial

    num = a**c / factorial(c)
    den = sum(a**k / factorial(k) for k in range(c + 1))
    return num / den


def erlang_c(c: int, a: float) -> float:
    """Probabilidad de esperar C(c, a) (Erlang C), requiere a < c.

    C(c, a) = B / (1 - ρ(1 - B)),  con B = B(c, a) y ρ = a/c.
    """
    if a >= c:
        raise ValueError("Erlang C exige a < c (sistema estable)")
    b = erlang_b(c, a)
    rho = a / c
    return b / (1.0 - rho * (1.0 - b))


def mmc(lam: float, mu: float, c: int) -> dict:
    """Medidas teóricas de M/M/c con espera (cola ilimitada, FIFO).

    Fórmulas:
        a   = λ/μ,  ρ = a/c
        C   = Erlang C(c, a)                 (= P(esperar), por PASTA)
        Lq  = C·ρ/(1-ρ)                      (pacientes en cola)
        Wq  = Lq/λ                           (días; Little aplicado a la cola)
        W   = Wq + 1/μ                       (días)
        L   = Lq + a                         (pacientes en el sistema)
    """
    a = carga_ofrecida(lam, mu)
    rho = a / c
    C = erlang_c(c, a)
    Lq = C * rho / (1.0 - rho)
    Wq = Lq / lam
    W = Wq + 1.0 / mu
    L = Lq + a
    return {
        "a": a, "rho": rho, "p_esperar": C, "Lq": Lq, "Wq": Wq, "W": W, "L": L,
        "utilizacion": rho, "p_bloqueo": 0.0, "lam_efectiva": lam,
        "Wq_h": Wq * 24.0, "W_h": W * 24.0,
    }


def mmcc(lam: float, mu: float, c: int) -> dict:
    """Medidas teóricas de M/M/c/c sin espera (pérdida pura, fórmula de Erlang B).

    Fórmulas:
        B          = B(c, a)  (recurrencia)             (= P(bloqueo), por PASTA)
        λ_ef       = λ(1 - B)                            (pacientes admitidos/día)
        L          = λ_ef/μ = a(1 - B)                   (camas ocupadas en promedio)
        utilización por cama = L/c
        W          = 1/μ (el que entra no espera), Lq = Wq = 0
    """
    a = carga_ofrecida(lam, mu)
    B = erlang_b(c, a)
    lam_ef = lam * (1.0 - B)
    L = lam_ef / mu
    return {
        "a": a, "p_bloqueo": B, "lam_efectiva": lam_ef, "L": L, "utilizacion": L / c,
        "W": 1.0 / mu, "Wq": 0.0, "Lq": 0.0, "p_esperar": 0.0,
        "W_h": 24.0 / mu, "Wq_h": 0.0,
    }


def little(L: float, lam: float, W: float) -> dict:
    """Ley de Little L = λ·W con λ = tasa de pacientes que *realmente entran* al sistema.

    Devuelve L, λ·W y el error relativo |L - λW| / L.
    """
    lw = lam * W
    err = abs(L - lw) / L if L != 0 else float("nan")
    return {"L": L, "lam_W": lw, "error_relativo": err, "cumple": isclose(L, lw, rel_tol=1e-9)}

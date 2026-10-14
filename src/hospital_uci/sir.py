"""Modelo SIR (Task 4.1a) y tasa de llegadas adicionales de pacientes críticos λ_epi(t).

Unidades: tiempo en DÍAS; β, γ en 1/día; S, I, R en personas; λ_epi y λ(t) en pacientes/día.

Ecuaciones (hoja, Task 4.1):
    dS/dt = -β S I / N
    dI/dt =  β S I / N - γ I
    dR/dt =  γ I
    λ_epi(t) = p_h · β · S(t) · I(t) / N          (p_h·(nuevos infectados por día))
    λ(t)     = λ + λ_epi(t)

Método: Runge–Kutta clásico de 4.º orden (RK4) con paso fijo. Es de un solo paso, explícito, de error global O(dt⁴),
y el sistema es suave y no rígido (β, γ ~ 0.3 y 0.2 por día), por lo que dt = 0.1 día ya da error ≪ 1e-6 relativo.
``integrar_referencia`` usa DOP853 con tolerancias estrictas solo como contraste independiente.

NO implementado aquí: generación de llegadas por adelgazamiento y simulación del hospital con λ(t) (Task 4.1b–c).
"""
from __future__ import annotations

import numpy as np

# Parámetros del grupo 5 (hoja, Task 4.1)
N_DEFECTO = 500_000
GAMMA_DEFECTO = 0.2
I0_DEFECTO = 100
R0_DEFECTO = 0
BETA_GRUPO5 = 0.30
P_H_DEFECTO = 0.0002
LAMBDA_BASE_GRUPO5 = 3.8


def _derivadas(y, beta, gamma, N):
    S, I, _ = y
    inf = beta * S * I / N
    return np.array([-inf, inf - gamma * I, gamma * I])


def lambda_epi(S, I, beta=BETA_GRUPO5, p_h=P_H_DEFECTO, N=N_DEFECTO):
    """λ_epi = p_h · β · S · I / N  (pacientes críticos adicionales por día). Acepta escalares o arreglos."""
    return p_h * beta * np.asarray(S) * np.asarray(I) / N


def integrar_sir(beta=BETA_GRUPO5, gamma=GAMMA_DEFECTO, N=N_DEFECTO, I0=I0_DEFECTO, R0=R0_DEFECTO,
                 dias=200.0, dt=0.1, p_h=P_H_DEFECTO) -> dict:
    """Integra el SIR con RK4 de paso fijo ``dt`` desde t=0 hasta t=dias.

    Devuelve dict con arreglos ``t, S, I, R, lambda_epi`` (longitud dias/dt + 1). Requiere dias/dt entero.
    """
    n = int(round(dias / dt))
    if abs(n * dt - dias) > 1e-9:
        raise ValueError("dias debe ser múltiplo de dt")
    y = np.array([N - I0 - R0, I0, R0], dtype=float)
    Y = np.empty((n + 1, 3))
    Y[0] = y
    for k in range(n):
        k1 = _derivadas(y, beta, gamma, N)
        k2 = _derivadas(y + 0.5 * dt * k1, beta, gamma, N)
        k3 = _derivadas(y + 0.5 * dt * k2, beta, gamma, N)
        k4 = _derivadas(y + dt * k3, beta, gamma, N)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        Y[k + 1] = y
    t = np.linspace(0.0, dias, n + 1)
    return {"t": t, "S": Y[:, 0], "I": Y[:, 1], "R": Y[:, 2],
            "lambda_epi": lambda_epi(Y[:, 0], Y[:, 1], beta, p_h, N)}


def integrar_referencia(beta=BETA_GRUPO5, gamma=GAMMA_DEFECTO, N=N_DEFECTO, I0=I0_DEFECTO, R0=R0_DEFECTO,
                        dias=200.0, dt=0.01, p_h=P_H_DEFECTO) -> dict:
    """Solución de referencia con DOP853 (rtol=1e-11, atol=1e-9), evaluada en la malla de paso ``dt``."""
    from scipy.integrate import solve_ivp

    n = int(round(dias / dt))
    t = np.linspace(0.0, dias, n + 1)
    sol = solve_ivp(lambda _t, y: _derivadas(y, beta, gamma, N), (0.0, dias), [N - I0 - R0, I0, R0],
                    method="DOP853", t_eval=t, rtol=1e-11, atol=1e-9)
    S, I, R = sol.y
    return {"t": t, "S": S, "I": I, "R": R, "lambda_epi": lambda_epi(S, I, beta, p_h, N)}


def pico(curva: dict, clave="lambda_epi") -> dict:
    """Día y valor del máximo de ``curva[clave]`` en la malla de la curva."""
    k = int(np.argmax(curva[clave]))
    return {"dia": float(curva["t"][k]), "valor": float(curva[clave][k])}


def exportar_csv(curva: dict, ruta: str, lam_base: float = LAMBDA_BASE_GRUPO5, cada: int = 1) -> None:
    """Escribe t,S,I,R,lambda_epi,lambda_total (pacientes/día) en CSV; ``cada`` submuestrea la malla."""
    import csv

    with open(ruta, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_dias", "S", "I", "R", "lambda_epi_pac_dia", "lambda_total_pac_dia"])
        for k in range(0, len(curva["t"]), cada):
            le = float(curva["lambda_epi"][k])
            w.writerow([f"{curva['t'][k]:.6f}", f"{curva['S'][k]:.6f}", f"{curva['I'][k]:.6f}",
                        f"{curva['R'][k]:.6f}", f"{le:.10f}", f"{lam_base + le:.10f}"])


def cargar_csv(ruta: str) -> dict:
    """Lee un CSV escrito por ``exportar_csv`` y devuelve arreglos numpy por columna."""
    datos = np.genfromtxt(ruta, delimiter=",", names=True)
    return {n: datos[n] for n in datos.dtype.names}


def lambda_epi_en(t, curva: dict):
    """Interpola linealmente λ_epi en instantes ``t`` (días) a partir de una curva/CSV cargado.

    Con la malla fina (dt = 0.1 día) del CSV ``sir_lambda_epi_dt0.1.csv`` el error de interpolación es despreciable.
    """
    clave = "lambda_epi" if "lambda_epi" in curva else "lambda_epi_pac_dia"
    tt = curva["t"] if "t" in curva else curva["t_dias"]
    return np.interp(t, tt, curva[clave])

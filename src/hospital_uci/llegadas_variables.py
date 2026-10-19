# Variante sin espera (M/Mt/c/c) con llegadas no homogéneas λ(t) por adelgazamiento y cierres temporales de camas.
# Se usa en los Tasks 3.1c y 4.1b-c; no modifica simulador.py (que sigue siendo el simulador validado del Task 1).
from __future__ import annotations

import heapq
import math

import numpy as np

from .simulador import ic_t, semillas_independientes


# Devuelve una función λ(t) constante (pacientes/día); sirve para validar el adelgazamiento contra Erlang B.
def tasa_constante(lam):
    return lambda t: lam


# Construye λ(t) = λ + λepi(t - t_inicio) a partir de la curva SIR del 4.1a (CSV con malla de 0.1 días).
# Antes de t_inicio la epidemia no existe (λepi = 0); dentro de la malla se interpola linealmente.
def tasa_ola(curva, lam_base=3.8, t_inicio=0.0):
    tt = np.asarray(curva["t_dias"] if "t_dias" in curva else curva["t"], dtype=float)
    le = np.asarray(curva["lambda_epi_pac_dia"] if "lambda_epi_pac_dia" in curva else curva["lambda_epi"], dtype=float)
    t_fin = float(tt[-1])

    def tasa(t):
        s = t - t_inicio
        if s < 0.0 or s > t_fin:
            return lam_base
        return lam_base + float(np.interp(s, tt, le))

    return tasa


# Estancia lognormal con media 1/μ y coeficiente de variación cs (fórmulas del Task 2.2a):
# s² = ln(1 + CS²) y m = ln(1/μ) − s²/2, donde m y s son la media y la desviación de ln(estancia).
def estancia_lognormal(mu, cs):
    s2 = math.log(1.0 + cs * cs)
    m = math.log(1.0 / mu) - s2 / 2.0
    s = math.sqrt(s2)
    return lambda rng: float(rng.lognormal(m, s))


# Capacidad vigente en t: c menos las camas cerradas por cada cierre (t0, t1, k_off) con t0 <= t < t1.
def capacidad_en(t, c, cierres):
    k = 0
    for t0, t1, k_off in cierres:
        if t0 <= t < t1:
            k += k_off
    return c - k


# Método de adelgazamiento (Lewis y Shedler): candidatos de un Poisson homogéneo de tasa λmax,
# y cada candidato en t se acepta con probabilidad λ(t)/λmax. Exige λ(t) <= λmax en todo t.
def llegadas_adelgazamiento(rng, tasa, lam_max, t_ini, t_fin):
    t = t_ini
    out = []
    while True:
        t += rng.exponential(1.0 / lam_max)
        if t > t_fin:
            break
        lt = tasa(t)
        if lt > lam_max * (1.0 + 1e-12):
            raise ValueError(f"λ({t:.3f}) = {lt:.4f} supera λmax = {lam_max:.4f}")
        if rng.random() * lam_max <= lt:
            out.append(t)
    return np.array(out)


# Corre UNA réplica de la variante sin espera con capacidad que puede bajar durante cierres.
# Llegadas: por adelgazamiento con (tasa, lam_max) o un arreglo fijo `llegadas` (para casos a mano).
# Al cerrar camas no se expulsa a nadie: solo se bloquea a quien llega con ocupadas >= capacidad(t).
# Empates: la salida se procesa antes que la llegada (misma convención que simulador.py).
def simular_perdida(c, horizonte, semilla=None, tasa=None, lam_max=None, mu=0.25, servicio=None,
                    cierres=(), llegadas=None, t_ini=0.0):
    if not (c >= 1 and int(c) == c):
        raise ValueError("c debe ser un entero >= 1")
    c = int(c)
    rng = semilla if isinstance(semilla, np.random.Generator) else np.random.default_rng(semilla)
    if llegadas is None:
        if tasa is None or lam_max is None:
            raise ValueError("se requiere tasa y lam_max, o un arreglo de llegadas")
        llegadas = llegadas_adelgazamiento(rng, tasa, lam_max, t_ini, horizonte)
    llegadas = np.asarray(llegadas, dtype=float)
    if servicio is None:
        escala = 1.0 / mu
        servicio = lambda g: g.exponential(escala)  # noqa: E731

    salidas = []                 # heap con los instantes de salida de los pacientes en cama
    cambio_t = [t_ini]           # instantes en que cambia la ocupación
    cambio_o = [0]               # ocupación vigente desde cada instante
    bloqueado = np.zeros(len(llegadas), dtype=bool)
    for i, t in enumerate(llegadas):
        while salidas and salidas[0] <= t:       # salida antes que llegada en empates
            ts = heapq.heappop(salidas)
            cambio_t.append(ts)
            cambio_o.append(cambio_o[-1] - 1)
        if cambio_o[-1] < capacidad_en(t, c, cierres):
            heapq.heappush(salidas, t + servicio(rng))
            cambio_t.append(t)
            cambio_o.append(cambio_o[-1] + 1)
        else:
            bloqueado[i] = True
    while salidas and salidas[0] <= horizonte:
        ts = heapq.heappop(salidas)
        cambio_t.append(ts)
        cambio_o.append(cambio_o[-1] - 1)
    return {"c": c, "cierres": tuple(cierres), "t_ini": t_ini, "horizonte": horizonte,
            "t_llegada": llegadas, "bloqueado": bloqueado,
            "cambio_t": np.array(cambio_t), "cambio_o": np.array(cambio_o)}


# Llegadas y bloqueados (cohorte por instante de llegada) en cada ventana [bordes[k], bordes[k+1]).
def conteo_por_ventana(res, bordes):
    bordes = np.asarray(bordes, dtype=float)
    k = np.searchsorted(bordes, res["t_llegada"], side="right") - 1
    ok = (k >= 0) & (k < len(bordes) - 1)
    n = np.bincount(k[ok], minlength=len(bordes) - 1)
    b = np.bincount(k[ok], weights=res["bloqueado"][ok].astype(float), minlength=len(bordes) - 1)
    return n.astype(int), b.astype(int)


# Fracción del tiempo de cada ventana con todas las camas vigentes ocupadas (ocupadas >= capacidad(t)).
# Por PASTA, una llegada de Poisson en t ve el estado de t, así que esto estima P(bloqueo) en esa ventana;
# también devuelve las camas ocupadas promedio por ventana (promedio ponderado por tiempo).
def tiempo_lleno_por_ventana(res, bordes):
    bordes = np.asarray(bordes, dtype=float)
    c, cierres = res["c"], res["cierres"]
    # puntos de quiebre: cambios de ocupación, extremos de cierres y bordes de ventanas
    extra = [x for t0, t1, _ in cierres for x in (t0, t1)]
    pts = np.unique(np.concatenate([res["cambio_t"], np.asarray(extra, dtype=float), bordes]))
    pts = pts[(pts >= res["t_ini"]) & (pts <= res["horizonte"])]
    # ocupación vigente en cada tramo [pts[j], pts[j+1]) (último cambio con instante <= pts[j])
    idx = np.searchsorted(res["cambio_t"], pts[:-1], side="right") - 1
    ocup = res["cambio_o"][idx]
    medio = 0.5 * (pts[:-1] + pts[1:])
    cap = np.array([capacidad_en(x, c, cierres) for x in medio])
    dur = np.diff(pts)
    F_lleno = np.concatenate([[0.0], np.cumsum(dur * (ocup >= cap))])
    F_ocup = np.concatenate([[0.0], np.cumsum(dur * ocup)])
    # integrales acumuladas evaluadas en los bordes (son lineales a trozos, así que interp es exacto)
    L = np.diff(bordes)
    lleno = np.diff(np.interp(bordes, pts, F_lleno)) / L
    ocupadas = np.diff(np.interp(bordes, pts, F_ocup)) / L
    return lleno, ocupadas


# Corre n_rep réplicas independientes (SeedSequence.spawn) y resume cada una con `resumen(res)`,
# para no guardar en memoria todas las llegadas: una observación por réplica para los IC.
def replicar_perdida(n_rep, semilla_maestra, resumen, **kwargs):
    hijas, estados = semillas_independientes(semilla_maestra, n_rep)
    out = []
    for h, e in zip(hijas, estados):
        r = resumen(simular_perdida(semilla=h, **kwargs))
        r["semilla_estado"] = e
        out.append(r)
    return out


# Media e IC t de 95 % entre réplicas, por columna, de una matriz (réplicas x ventanas).
def ic_por_columna(matriz, nivel=0.95):
    m = np.asarray(matriz, dtype=float)
    filas = [ic_t(m[:, j], nivel) for j in range(m.shape[1])]
    return {k: np.array([f[k] for f in filas]) for k in ("media", "ic_inf", "ic_sup", "semiancho")}

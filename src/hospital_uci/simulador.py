"""Simulador de eventos discretos CORREGIDO y validado: M/M/c (con espera) y M/M/c/c (sin espera).

Reemplaza a ``anexo/B_simulador_original_ia.py`` (que NO debe importarse). La auditoría contra los siete
puntos de la hoja (Task 1.3) está en ``notebooks/task1_teoria_simulador_y_validacion.ipynb``.

Unidades: tiempos en DÍAS, tasas en 1/día. Wq y W se devuelven en días (multiplique por 24 para horas).

Diseño
------
* Lista de eventos futuros (FEL) = heap ``(tiempo, prioridad, secuencia, tipo, paciente)``.
  Empates: la salida se procesa antes que la llegada (prioridad 0 < 1), de modo que una cama que se libera en el
  mismo instante de una llegada se considera libre.
* Las llegadas se generan hasta ``horizonte``. Después de eso solo quedan salidas ("drenado"): los pacientes de
  la cohorte (los que llegaron en [calentamiento, horizonte]) se siguen hasta que salen, de modo que Wq y W
  no se sesgan por censura (punto 6 de la hoja).
* L, Lq y camas ocupadas son promedios PONDERADOS POR TIEMPO sobre la ventana [calentamiento, horizonte]
  (punto 1). Los pacientes que están dentro al llegar el horizonte cuentan en esos promedios hasta ese instante.
* Wq se mide desde la llegada del paciente hasta que empieza su servicio (punto 7); W desde la llegada hasta la salida.
* p_esperar y p_bloqueo = fracción de llegadas de la cohorte que encontró todas las camas ocupadas
  (espera / bloqueo según la variante). Por PASTA coinciden con Erlang C / Erlang B.
* Réplicas: ``replicar`` deriva semillas con ``numpy.random.SeedSequence.spawn`` (flujos independientes, punto 4)
  y ``resumir`` calcula el IC t de 95 % usando UNA observación por réplica (punto 2).

Ganchos de extensión YA implementados (y probados solo con casos deterministas y exponenciales):
    ``interarribo(rng) -> float`` y ``servicio(rng) -> float`` reemplazan a las exponenciales por defecto.
NO implementados: cierre parcial de camas (capacidad variable), llegadas no homogéneas λ(t), distribuciones
lognormal/determinista predefinidas, series temporales de bloqueo. Ver la sección "Extensiones necesarias" del ``README.md``.
"""
from __future__ import annotations

import heapq
import math
from collections import deque
from typing import Callable, Optional, Sequence

import numpy as np
from scipy import stats

_SALIDA, _LLEGADA = 0, 1  # prioridad de desempate: salida antes que llegada


def simular(
    lam: float,
    mu: float,
    c: int,
    horizonte: float,
    calentamiento: float = 0.0,
    semilla=None,
    con_espera: bool = True,
    muestreo_dt: Optional[float] = None,
    interarribo: Optional[Callable[[np.random.Generator], float]] = None,
    servicio: Optional[Callable[[np.random.Generator], float]] = None,
) -> dict:
    """Corre UNA réplica y devuelve un diccionario de medidas.

    Entradas
        lam            λ, pacientes/día (solo se usa si ``interarribo`` es None)
        mu             μ, 1/día por cama (solo se usa si ``servicio`` es None); estancia media 1/μ
        c              camas (entero >= 1)
        horizonte      instante final de llegadas y de la ventana de promedios (días)
        calentamiento  instante desde el cual se recolectan estadísticas (días); 0 <= calentamiento < horizonte
        semilla        int, ``SeedSequence`` o ``Generator`` de numpy (None = aleatoria)
        con_espera     True: M/M/c (cola FIFO ilimitada). False: M/M/c/c (el que no encuentra cama se pierde)
        muestreo_dt    si se da, devuelve además ``traza`` con la cola y el sistema en t = 0, dt, 2dt, ... <= horizonte
        interarribo/servicio  ganchos opcionales (ver docstring del módulo)

    Salidas (claves)
        L, Lq            pacientes en sistema / en cola, promedio ponderado por tiempo en [calentamiento, horizonte]
        camas_ocupadas   promedio ponderado de camas ocupadas;  utilizacion = camas_ocupadas / c
        Wq, W            espera en cola y tiempo total en el sistema (DÍAS), promedio sobre los pacientes ADMITIDOS de la cohorte
        p_esperar, p_bloqueo   fracción de llegadas de la cohorte que esperan / son bloqueadas
        lam_efectiva     admitidos de la cohorte / (horizonte - calentamiento), pacientes/día
        n_llegadas, n_admitidos, n_bloqueados  conteos de la cohorte
        n_dentro_al_horizonte  pacientes en el sistema (servicio + cola) en el instante ``horizonte``
        semilla_entropia  entropía (int) del generador, útil para reportar semillas
    """
    if not (c >= 1 and int(c) == c):
        raise ValueError("c debe ser un entero >= 1")
    if not (0 <= calentamiento < horizonte):
        raise ValueError("se requiere 0 <= calentamiento < horizonte")
    if interarribo is None and lam <= 0:
        raise ValueError("lam debe ser > 0")
    if servicio is None and mu <= 0:
        raise ValueError("mu debe ser > 0")
    c = int(c)

    rng = semilla if isinstance(semilla, np.random.Generator) else np.random.default_rng(semilla)
    if interarribo is None:
        escala_ll = 1.0 / lam
        interarribo = lambda g: g.exponential(escala_ll)  # noqa: E731
    if servicio is None:
        escala_sv = 1.0 / mu
        servicio = lambda g: g.exponential(escala_sv)  # noqa: E731

    fel: list = []
    seq = 0

    def programar(t, prio, tipo, pid):
        nonlocal seq
        heapq.heappush(fel, (t, prio, seq, tipo, pid))
        seq += 1

    t_ll = interarribo(rng)
    siguiente_id = 0
    if t_ll <= horizonte:
        programar(t_ll, _LLEGADA, _LLEGADA, siguiente_id)
        siguiente_id += 1

    ocupadas = 0
    cola: deque = deque()
    t_llegada: dict = {}      # id -> instante de llegada (solo pacientes que siguen en el sistema)
    en_cohorte: set = set()   # ids de pacientes admitidos cuya llegada fue >= calentamiento

    ultimo = 0.0
    area_L = area_Lq = area_ocup = 0.0
    n_llegadas = n_esperaron = n_bloqueados = n_admitidos = 0
    sum_wq = sum_w = 0.0
    n_servidos_cohorte = 0
    n_dentro_horizonte = None

    traza_t = traza_q = traza_n = None
    sig_muestra = 0.0
    if muestreo_dt is not None:
        traza_t, traza_q, traza_n = [], [], []

    def acumular(t_hasta):
        """Integra el estado vigente (constante) sobre [ultimo, t_hasta] ∩ [calentamiento, horizonte]."""
        nonlocal area_L, area_Lq, area_ocup
        a0 = max(ultimo, calentamiento)
        a1 = min(t_hasta, horizonte)
        if a1 > a0:
            dt = a1 - a0
            area_L += (ocupadas + len(cola)) * dt
            area_Lq += len(cola) * dt
            area_ocup += ocupadas * dt

    def muestrear_hasta(t_evento):
        """Registra el estado vigente en los instantes de la malla estrictamente anteriores al evento."""
        nonlocal sig_muestra
        while sig_muestra < t_evento and sig_muestra <= horizonte:
            traza_t.append(sig_muestra)
            traza_q.append(len(cola))
            traza_n.append(ocupadas + len(cola))
            sig_muestra += muestreo_dt

    horizonte_cerrado = False
    while fel:
        t, _, _, tipo, pid = heapq.heappop(fel)

        # al cruzar el horizonte: cerrar el último tramo y fotografiar el estado (aún sin aplicar este evento)
        if t > horizonte and not horizonte_cerrado:
            acumular(horizonte)
            if muestreo_dt is not None:
                muestrear_hasta(math.inf)
            n_dentro_horizonte = ocupadas + len(cola)
            horizonte_cerrado = True
        if not horizonte_cerrado:
            acumular(t)
            if muestreo_dt is not None:
                muestrear_hasta(t)
        ultimo = t

        if tipo == _LLEGADA:
            sig = t + interarribo(rng)
            if sig <= horizonte:
                programar(sig, _LLEGADA, _LLEGADA, siguiente_id)
                siguiente_id += 1
            cuenta = t >= calentamiento
            if cuenta:
                n_llegadas += 1
            if ocupadas < c:
                ocupadas += 1
                t_llegada[pid] = t
                if cuenta:
                    en_cohorte.add(pid)
                    n_admitidos += 1
                    # espera 0: empieza servicio al llegar
                programar(t + servicio(rng), _SALIDA, _SALIDA, pid)
            elif con_espera:
                cola.append(pid)
                t_llegada[pid] = t
                if cuenta:
                    en_cohorte.add(pid)
                    n_admitidos += 1
                    n_esperaron += 1
            else:
                if cuenta:
                    n_bloqueados += 1
        else:  # salida
            if pid in en_cohorte:
                sum_w += t - t_llegada[pid]
                n_servidos_cohorte += 1
                en_cohorte.discard(pid)
            del t_llegada[pid]
            if cola:
                nuevo = cola.popleft()
                if nuevo in en_cohorte:
                    sum_wq += t - t_llegada[nuevo]   # espera medida DESDE LA LLEGADA (punto 7)
                programar(t + servicio(rng), _SALIDA, _SALIDA, nuevo)
            else:
                ocupadas -= 1
        # nota: en el cálculo de sum_wq solo entran admitidos de la cohorte; los de espera 0 aportan 0.

    if not horizonte_cerrado:  # nunca hubo eventos posteriores al horizonte
        acumular(horizonte)
        if muestreo_dt is not None:
            muestrear_hasta(math.inf)
        n_dentro_horizonte = ocupadas + len(cola)

    T = horizonte - calentamiento
    # Todos los admitidos de la cohorte ya salieron (drenado), así que n_servidos_cohorte == n_admitidos.
    res = {
        "L": area_L / T,
        "Lq": area_Lq / T,
        "camas_ocupadas": area_ocup / T,
        "utilizacion": area_ocup / (c * T),
        "Wq": sum_wq / n_admitidos if n_admitidos else float("nan"),
        "W": sum_w / n_servidos_cohorte if n_servidos_cohorte else float("nan"),
        "p_esperar": n_esperaron / n_llegadas if n_llegadas else float("nan"),
        "p_bloqueo": n_bloqueados / n_llegadas if n_llegadas else float("nan"),
        "lam_efectiva": n_admitidos / T,
        "n_llegadas": n_llegadas,
        "n_admitidos": n_admitidos,
        "n_bloqueados": n_bloqueados,
        "n_dentro_al_horizonte": n_dentro_horizonte,
    }
    if muestreo_dt is not None:
        res["traza"] = {"t": np.array(traza_t), "cola": np.array(traza_q), "sistema": np.array(traza_n)}
    return res


def semillas_independientes(semilla_maestra: int, n: int):
    """Devuelve ``n`` SeedSequence hijas independientes (``SeedSequence.spawn``) y su entropía representativa."""
    hijas = np.random.SeedSequence(semilla_maestra).spawn(n)
    return hijas, [int(h.generate_state(1)[0]) for h in hijas]


def replicar(n_rep: int, semilla_maestra: int, **kwargs) -> list:
    """Corre ``n_rep`` réplicas independientes de ``simular`` (kwargs: lam, mu, c, horizonte, ...).

    Cada réplica usa una SeedSequence hija distinta de ``semilla_maestra`` (punto 4 de la hoja).
    Devuelve la lista de diccionarios; cada uno incluye ``semilla_estado`` (primer entero de su estado).
    """
    hijas, estados = semillas_independientes(semilla_maestra, n_rep)
    out = []
    for h, e in zip(hijas, estados):
        r = simular(semilla=h, **kwargs)
        r["semilla_estado"] = e
        out.append(r)
    return out


def ic_t(valores: Sequence[float], nivel: float = 0.95) -> dict:
    """IC t de ``nivel`` para la media usando UNA observación por réplica (punto 2): n-1 g.l."""
    x = np.asarray(valores, dtype=float)
    x = x[~np.isnan(x)]
    n = len(x)
    m = float(x.mean())
    if n < 2:
        return {"media": m, "ic_inf": float("nan"), "ic_sup": float("nan"), "semiancho": float("nan"), "n": n, "desv": float("nan")}
    s = float(x.std(ddof=1))
    h = float(stats.t.ppf(0.5 + nivel / 2.0, n - 1) * s / math.sqrt(n))
    return {"media": m, "ic_inf": m - h, "ic_sup": m + h, "semiancho": h, "n": n, "desv": s}


def resumir(replicas: list, claves=("L", "Lq", "camas_ocupadas", "utilizacion", "Wq", "W", "p_esperar", "p_bloqueo", "lam_efectiva")) -> dict:
    """Media e IC 95 % entre réplicas para cada clave (una observación por réplica)."""
    return {k: ic_t([r[k] for r in replicas]) for k in claves}

"""
Simulador de eventos discretos del hospital: M/M/c (con espera) y M/M/c/c (sin espera).

Primera version generada por IA (Claude Sonnet 5.5) en Claude Code.
Usa una lista de eventos futuros (heap). Tiempos en dias.
"""

import heapq
import random
import statistics
import math

from scipy import stats

ARRIVAL = 0
DEPARTURE = 1


def simulate(lam, mu, c, horizon, warmup=0.0, seed=None, wait=True):
    """Corre una replica.

    lam: tasa de llegadas (pacientes/dia)
    mu: tasa de servicio por cama (1/dia)
    c: numero de camas
    horizon: tiempo total simulado (dias)
    warmup: tiempo inicial que se descarta (dias)
    seed: semilla
    wait: True -> M/M/c (con espera); False -> M/M/c/c (sin espera, se bloquea)
    """
    rng = random.Random(seed)

    fel = []  # lista de eventos futuros: (tiempo, tipo, id)
    heapq.heappush(fel, (rng.expovariate(lam), ARRIVAL, 0))
    next_id = 1

    busy = 0
    queue = []  # FIFO de (tiempo_llegada, cuenta_para_estadisticas)

    last = 0.0
    area_L = 0.0
    area_Lq = 0.0
    area_busy = 0.0

    n_arrivals = 0
    n_waited = 0
    n_blocked = 0
    waits = []

    while fel:
        t, kind, pid = heapq.heappop(fel)

        if t > horizon:
            # cerrar el ultimo tramo hasta el horizonte
            if horizon > warmup:
                dt = horizon - max(last, warmup)
                if dt > 0:
                    area_L += (busy + len(queue)) * dt
                    area_Lq += len(queue) * dt
                    area_busy += busy * dt
            break

        # acumular areas del tramo anterior (solo despues del calentamiento)
        if t > warmup:
            dt = t - max(last, warmup)
            area_L += (busy + len(queue)) * dt
            area_Lq += len(queue) * dt
            area_busy += busy * dt
        last = t

        counting = t >= warmup

        if kind == ARRIVAL:
            # programar la siguiente llegada
            heapq.heappush(fel, (t + rng.expovariate(lam), ARRIVAL, next_id))
            next_id += 1

            if counting:
                n_arrivals += 1

            if busy < c:
                busy += 1
                heapq.heappush(fel, (t + rng.expovariate(mu), DEPARTURE, pid))
                if counting:
                    waits.append(0.0)
            elif wait:
                queue.append((t, counting))
                if counting:
                    n_waited += 1
            else:
                if counting:
                    n_blocked += 1
        else:
            if queue:
                arr_t, was_counting = queue.pop(0)
                heapq.heappush(fel, (t + rng.expovariate(mu), DEPARTURE, pid))
                if was_counting:
                    waits.append(t - arr_t)
            else:
                busy -= 1

    T = horizon - warmup
    return {
        "L": area_L / T,
        "Lq": area_Lq / T,
        "Wq": statistics.mean(waits) if waits else 0.0,
        "utilizacion": area_busy / (c * T),
        "p_esperar": n_waited / n_arrivals if n_arrivals else 0.0,
        "p_bloqueo": n_blocked / n_arrivals if n_arrivals else 0.0,
    }


def replicate(lam, mu, c, horizon, warmup, n_reps=30, base_seed=1000, wait=True):
    """Corre n_reps replicas y devuelve media e IC 95% por medida."""
    results = []
    for i in range(n_reps):
        results.append(simulate(lam, mu, c, horizon, warmup, seed=base_seed + i, wait=wait))

    summary = {}
    for key in results[0]:
        vals = [r[key] for r in results]
        m = statistics.mean(vals)
        s = statistics.stdev(vals)
        h = stats.t.ppf(0.975, n_reps - 1) * s / math.sqrt(n_reps)
        summary[key] = (m, m - h, m + h)
    return summary


if __name__ == "__main__":
    lam, mu, c = 3.8, 0.25, 18
    print("Con espera:", replicate(lam, mu, c, 2000, 200, wait=True))
    print("Sin espera:", replicate(lam, mu, c, 2000, 200, wait=False))

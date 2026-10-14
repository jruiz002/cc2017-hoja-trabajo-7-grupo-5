# Hoja de Trabajo 7 — Grupo 5 (CC2017 Modelación y Simulación, UVG)

Hospital Regional del Altiplano: la UCI se modela como una cola M/M/c (con espera) y M/M/c/c (sin espera), con un simulador de eventos discretos validado contra la teoría de colas y un modelo de epidemia SIR acoplado en un solo sentido.

Parámetros del grupo 5: **λ = 3.8 pacientes/día, μ = 0.25/día, c = 18 camas, a = 15.2, ρ = 0.8444, k_off = 6 camas (durante 7 días), β = 0.30/día.**

## Estado por Task

| Task | Contenido | Notebook | Estado |
|---|---|---|---|
| 1 (1.1–1.4) | Teoría, simulador generado por IA (anexo), auditoría y corrección, validación | `notebooks/task1_teoria_simulador_y_validacion.ipynb` | Hecho |
| 2 (2.1–2.4) | Estancia no exponencial, cierre parcial de camas, balance de predicciones | `notebooks/task2_estancias_y_cierre_de_camas.ipynb` | Pendiente (plantilla) |
| 3 (3.1–3.3) | Recomendación de camas, robustez y defensa | `notebooks/task3_recomendacion_de_camas.ipynb` | Pendiente (plantilla) |
| 4.1a | Modelo SIR, λepi(t) y pico de la ola | `notebooks/task4_epidemia_sir.ipynb` | Hecho |
| 4.1b–c | Simulación con λ(t) y camas de contingencia | `notebooks/task4_epidemia_sir.ipynb` | Pendiente |
| 4.1d | Límites del modelo de un solo sentido | `notebooks/task4_epidemia_sir.ipynb` | Pendiente |

## Cómo ejecutar

Requiere Python ≥ 3.11 (probado con 3.14).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest -q                      # 23 pruebas, ~2 s
.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/task4_epidemia_sir.ipynb                       # ~10 s
.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/task1_teoria_simulador_y_validacion.ipynb      # ~6 min
cd anexo && shasum -a 256 -c SHA256SUMS.txt        # el código original de la IA sigue intacto
```

También se pueden abrir los notebooks en VS Code o JupyterLab con el kernel de `.venv`. Cada notebook guarda sus resultados en `resultados/`.

## Estructura

| Ruta | Contenido |
|---|---|
| `src/hospital_uci/teoria.py` | Erlang B (recurrencia), Erlang C, M/M/c, M/M/c/c, ley de Little, parámetros del grupo |
| `src/hospital_uci/simulador.py` | **Simulador corregido y validado** (`simular`, `replicar`, `resumir`, `ic_t`) |
| `src/hospital_uci/sir.py` | SIR (RK4), λepi(t), pico, exportar/cargar CSV, interpolación |
| `notebooks/` | Un notebook por Task (ver tabla de estado) |
| `tests/` | Pruebas con respuesta conocida: casos D/D/c calculados a mano, M/M/1, Erlang B, SIR |
| `resultados/figuras/`, `resultados/datos/` | Figuras PNG, resultados en JSON, semillas y CSV del SIR |
| `anexo/` | Anexo del Task 1.2: prompt (`A_…`), código original de la IA **sin modificar** y de solo lectura (`B_…`), afirmaciones de la IA (`C_…`) y hashes |

**No se debe importar nada de `anexo/`**: ese simulador tiene un sesgo conocido (censura de Wq) y un desempate de eventos incorrecto; se conserva solo como evidencia.

## Uso mínimo del simulador

```python
import sys; sys.path.insert(0, "src")
from hospital_uci import simular, replicar, resumir, mmc, mmcc

# una réplica (tiempos en DÍAS; Wq y W se multiplican por 24 para horas)
r = simular(lam=3.8, mu=0.25, c=18, horizonte=11500, calentamiento=100, semilla=1, con_espera=False)

# 30 réplicas independientes + IC de 95 % entre réplicas
reps = replicar(30, 5004, lam=3.8, mu=0.25, c=18, horizonte=11500, calentamiento=100)
print(resumir(reps)["Lq"])     # {'media', 'ic_inf', 'ic_sup', 'semiancho', 'n', 'desv'}
```

Salidas por réplica: `L, Lq, camas_ocupadas, utilizacion, Wq, W, p_esperar, p_bloqueo, lam_efectiva, n_llegadas, n_admitidos, n_bloqueados, n_dentro_al_horizonte`.
Convenciones: promedios ponderados por tiempo en [calentamiento, horizonte]; Wq, W y probabilidades por cohorte de llegadas en esa ventana, seguida hasta que sale; en empates, la salida se procesa antes que la llegada.
Validación (Task 1.4): horizonte 11 500 días, calentamiento 100 días, 30 réplicas; las 12 medidas comparadas contienen el valor teórico en su IC de 95 %.

## Extensiones necesarias

Lo que **existe** está implementado y probado; lo que **falta** no está implementado todavía.

**Task 2**
- Existe: `teoria.mmc`, `mmcc`, `erlang_b` (para 2.1c: con estancia +25 %, μ' = μ/1.25 = 0.2 y a' = 19); `simular(..., servicio=...)`, un gancho que reemplaza la distribución de la estancia (probado solo con valores constantes y la exponencial); `muestreo_dt` (traza de la cola por réplica).
- Falta: generadores de estancia determinista y lognormal con CS = 2 (`s² = ln(1 + CS²)`, `m = ln(1/μ) − s²/2`; validar que la media simulada sea 1/μ); eventos de cierre y reapertura de camas (al bajar la capacidad, los pacientes en servicio no se expulsan); promedio de la cola entre réplicas en función del tiempo y tiempo de recuperación.
- Precaución: con estancias no exponenciales siguen valiendo la ley de Little y los IC, pero Erlang C ya no es exacto (Erlang B sí es insensible a la distribución de la estancia).

**Task 3**
- Existe: `simulador.simular` / `replicar` con `c` como parámetro (probado con c = 1 y c = 18); `teoria.erlang_b` y `mmc` para buscar el menor c con la teoría.
- Falta: el barrido de camas y su verificación por simulación; el cierre de camas del inciso 3.1c (misma extensión que el inciso 2.3; coordinar para no duplicar).

**Task 4.1b–c**
- Existe: `resultados/datos/sir_lambda_epi_dt0.1.csv` y `sir_lambda_epi_diario.csv` (columnas `t_dias, S, I, R, lambda_epi_pac_dia, lambda_total_pac_dia`), `sir.cargar_csv` y `sir.lambda_epi_en` (interpolación lineal). Pico de λepi: 1.3213 pac/día en el día 68.6; λ(t) máxima: 5.1213 pac/día (usar λmax ≥ 5.1213 en el adelgazamiento).
- Falta: llegadas no homogéneas por adelgazamiento (el gancho `interarribo(rng)` no recibe el reloj; hace falta pasarle el tiempo o añadir una tasa λ(t)); probabilidad de bloqueo por ventanas de tiempo y total de desviados con IC (100 réplicas); búsqueda de camas de contingencia para ≤ 5 % de desviados.
- Precaución: el calentamiento de 100 días se justificó para λ constante; con λ(t) hay que decidir qué hacer con el arranque con el hospital vacío.

Si se modifica `simulador.py`, añadir una prueba con respuesta conocida y volver a correr `pytest -q` antes de usar los resultados.

## Reproducibilidad

Todas las semillas están fijadas (maestras 5001–5005, 600000+j y 700000+j; ver el notebook de Task 1 y `resultados/datos/semillas_validacion.csv`). Con las mismas versiones de numpy y scipy se reproducen los números; con otras versiones, los valores de réplicas individuales podrían variar ligeramente, pero no las conclusiones.

## Autores

José Gerardo Ruiz García (23719), Gerardo André Fernández Cruz (23763) y Juan Diego Solís Martínez (23720).

# Hoja de Trabajo 7 — Grupo 5 (CC2017 Modelación y Simulación, UVG)

Hospital Regional del Altiplano: la UCI se modela como una cola M/M/c (con espera) y M/M/c/c (sin espera), con un simulador de eventos discretos validado contra la teoría de colas y un modelo de epidemia SIR acoplado en un solo sentido.

Parámetros del grupo 5: **λ = 3.8 pacientes/día, μ = 0.25/día, c = 18 camas, a = 15.2, ρ = 0.8444, k_off = 6 camas (durante 7 días), β = 0.30/día.**

## Estado por Task

| Task | Contenido | Notebook | Estado |
|---|---|---|---|
| 1 (1.1–1.4) | Teoría, simulador generado por IA (anexo), auditoría y corrección, validación | `notebooks/task1_teoria_simulador_y_validacion.ipynb` | Hecho |
| 2 (2.1–2.4) | Estancia no exponencial, cierre parcial de camas, balance de predicciones | `notebooks/task2_estancias_y_cierre_de_camas.ipynb` | Pendiente (plantilla) |
| 3 (3.1–3.3) | Recomendación de camas, robustez y defensa | `notebooks/task3_recomendacion_de_camas.ipynb` | Hecho (3.1 en el notebook; 3.2 y 3.3 en el informe) |
| 4.1a | Modelo SIR, λepi(t) y pico de la ola | `notebooks/task4_epidemia_sir.ipynb` | Hecho |
| 4.1b–c | Simulación con λ(t) y camas de contingencia | `notebooks/task4_epidemia_sir.ipynb` | Hecho |
| 4.1d | Límites del modelo de un solo sentido | `notebooks/task4_epidemia_sir.ipynb` | Pendiente |

## Cómo ejecutar

Requiere Python ≥ 3.11 (probado con 3.14).

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest -q                      # 32 pruebas, ~3 s
.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/task4_epidemia_sir.ipynb                       # ~10 s
.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/task1_teoria_simulador_y_validacion.ipynb      # ~6 min
.venv/bin/jupyter nbconvert --to notebook --execute --inplace notebooks/task3_recomendacion_de_camas.ipynb             # ~3 min (antes que el de Task 4)
cd anexo && shasum -a 256 -c SHA256SUMS.txt        # el código original de la IA sigue intacto
```

También se pueden abrir los notebooks en VS Code o JupyterLab con el kernel de `.venv`. Cada notebook guarda sus resultados en `resultados/`.

## Estructura

| Ruta | Contenido |
|---|---|
| `src/hospital_uci/teoria.py` | Erlang B (recurrencia), Erlang C, M/M/c, M/M/c/c, ley de Little, parámetros del grupo |
| `src/hospital_uci/simulador.py` | **Simulador corregido y validado** (`simular`, `replicar`, `resumir`, `ic_t`) |
| `src/hospital_uci/sir.py` | SIR (RK4), λepi(t), pico, exportar/cargar CSV, interpolación |
| `src/hospital_uci/llegadas_variables.py` | Variante sin espera con llegadas λ(t) por adelgazamiento, cierres o camas temporales por ventana de tiempo, P(bloqueo) por ventana y estancia lognormal (Tasks 3.1c y 4.1b–c) |
| `notebooks/` | Un notebook por Task (ver tabla de estado) |
| `tests/` | Pruebas con respuesta conocida: casos D/D/c calculados a mano, M/M/1, Erlang B, SIR, adelgazamiento y cierres |
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

**Task 3** (hecho)
- `simulador.simular` / `replicar` para verificar el barrido teórico de camas (3.1a y 3.1b); `llegadas_variables.estancia_lognormal` para CS = 2 (3.1b iii) con el gancho `servicio`; `llegadas_variables.simular_perdida(..., cierres=[(t0, t1, k_off)])` para el cierre del 3.1c (variante sin espera; sin expulsar pacientes).
- El cierre con cola (variante con espera) es del inciso 2.3 y no se implementó aquí.

**Task 4.1b–c** (hecho)
- `llegadas_variables.tasa_ola` y `llegadas_adelgazamiento` (λmax = 5.1213; se detiene si λ(t) > λmax); `tiempo_lleno_por_ventana` (P(bloqueo) por ventana, por PASTA) y `conteo_por_ventana` (desviados y llegadas). Camas de contingencia = "cierre" con `k_off` negativo durante la ola.
- Arranque: 100 días previos con λ = 3.8 antes del día 0 de la epidemia. "Durante la ola" = días con λepi ≥ 10 % de su pico (31.5 a 112.5).
- Validación (`tests/test_llegadas_variables.py`): casos a mano con cierre y camas temporales, Erlang B con λ constante por adelgazamiento y número de llegadas = integral de λ(t) escalonada.

Si se modifica `simulador.py`, añadir una prueba con respuesta conocida y volver a correr `pytest -q` antes de usar los resultados.

## Reproducibilidad

Todas las semillas están fijadas (maestras 5001–5005, 600000+j y 700000+j; ver el notebook de Task 1 y `resultados/datos/semillas_validacion.csv`). Task 3: 8101–8199 y 8301–8399 (`resultados/datos/task3_semillas.csv`); Task 4.1b–c: 8401 y 8402. Con las mismas versiones de numpy y scipy se reproducen los números; con otras versiones, los valores de réplicas individuales podrían variar ligeramente, pero no las conclusiones.

## Autores

José Gerardo Ruiz García (23719), Gerardo André Fernández Cruz (23763) y Juan Diego Solís Martínez (23720).

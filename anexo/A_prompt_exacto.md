# Anexo A — Prompt utilizado (Task 1.2)

## Parámetros y comprobaciones iniciales

Para el grupo g = 5, calcula a partir de las fórmulas de la hoja λ, c, a, ρ, k_off y β. Comprueba tus cálculos antes de continuar: deben resultar **λ = 3.8 pacientes/día, c = 18 camas, a = 15.2, ρ ≈ 0.8444, k_off = 6 camas y β = 0.30 por día**. Usa μ = 0.25 por día. Coloca los parámetros requeridos al inicio del texto del informe.

## Orden de trabajo: conserva primero el original

Trabaja por etapas. **Genera y guarda la primera versión del simulador antes de revisarla o corregirla.** Conserva esa versión literalmente, en un archivo distinto del código corregido. Conserva también una copia exacta de este prompt; indícame que debo guardar personalmente el texto enviado si tu interfaz no puede recuperarlo con fidelidad. Si afirmas que el código funciona, coincide con la teoría o superó pruebas, guarda esas afirmaciones textuales para contrastarlas después. No fabriques afirmaciones que no hayas hecho.

Una vez archivada la primera versión, continúa con la revisión, corrección y validación. No cambies el original retrospectivamente para que parezca mejor o peor. Explica qué verificaste realmente mediante ejecución y qué debo revisar y poder defender yo.

## Task 1.1: teoría de colas

Desarrolla las variantes **M/M/c con espera** y **M/M/c/c sin espera**. Muestra las fórmulas y evalúalas con código.

- Con espera: carga ofrecida a, utilización ρ, probabilidad de esperar de Erlang C, longitud media de cola Lq, tiempo medio de espera Wq y tiempo medio total W. Expresa los tiempos en horas.
- Sin espera: calcula la probabilidad de bloqueo de Erlang B mediante recurrencia, la tasa efectiva de ingresos y el número medio de camas ocupadas.
- Comprueba la ley de Little en ambas variantes. Para la variante sin espera usa la **tasa efectiva de ingresos** y explica en una línea el error que produciría usar λ.

Identifica claramente unidades, fórmulas, sustitución de valores y resultados numéricos.

## Task 1.2: primera versión generada por IA

Genera un simulador de eventos discretos para ambas variantes que use una **lista de eventos** y acepte como entradas λ, μ, c, horizonte, calentamiento y semilla. Debe reportar L, Lq, Wq, utilización, probabilidad de esperar y probabilidad de bloqueo. Guarda el código de esta primera versión sin modificaciones para el anexo. Después de guardarlo, puedes crear una versión corregida en otro archivo.

## Task 1.3: revisión y corrección

Audita el código original contra los **siete puntos exactos** de la hoja:

1. L y Lq se promedian ponderando cada estado por el tiempo que permaneció vigente.
2. Los intervalos de confianza usan la réplica como unidad de muestreo, no cada paciente.
3. Se descarta y justifica un periodo de calentamiento.
4. Las réplicas usan semillas independientes.
5. La variante sin espera usa la tasa efectiva de ingresos donde corresponde.
6. Los pacientes que siguen dentro al terminar el horizonte reciben un tratamiento consistente.
7. La espera se mide desde la llegada del paciente, no desde el inicio del servicio.

Presenta una tabla con una fila por punto y columnas para **cumplimiento, evidencia concreta, supuesto afectado si falla y prueba que lo detectaría**. Incluye por lo menos una prueba diseñada con una respuesta conocida de antemano. No inventes fallas: demuestra también los puntos que sí se cumplen. Entrega el código corregido por separado y redacta una evaluación honesta de la confiabilidad del original y de los errores que las pruebas permitieron descubrir.

## Task 1.4: validación reproducible

Ejecuta **al menos 30 réplicas independientes por variante** con el simulador corregido. Elige y justifica horizonte y calentamiento mediante una gráfica del **promedio móvil de la longitud de cola**. Para cada medida aplicable, reporta el valor teórico, la media simulada, el intervalo de confianza del 95 % calculado entre réplicas y si contiene el valor teórico.

Comprueba la ley de Little usando las medidas **obtenidas de la simulación**, independientemente de la comparación con la teoría, y reporta el error relativo. Contrasta tus afirmaciones iniciales sobre el simulador con las pruebas efectuadas. Si todas se confirman, diseña una prueba adicional que intente encontrar fallas y reporta lo que ocurrió. Muestra las semillas, el procedimiento y cualquier limitación; no inventes ejecuciones ni resultados.

## Task 4.1a: modelo SIR

Integra el modelo SIR durante **200 días** con N = 500 000, γ = 0.2 por día, I(0) = 100, R(0) = 0, S(0) = 499 900, β = 0.30 por día y pₕ = 0.0002. Usa las ecuaciones y la fórmula de λepi(t) indicadas en la hoja. Explica el método numérico elegido y comprueba que reducir el paso temporal no cambia materialmente el pico.

Grafica **λepi(t) junto con la tasa base λ = 3.8 pacientes/día**. Reporta el día y valor del pico de llegadas adicionales. Exporta los datos de la curva en un formato sencillo, como CSV, para que Juan Diego pueda utilizarlos más adelante. Este inciso debe quedar terminado sin esperar la recomendación de camas del Task 3. No simules todavía los incisos 4.1b–c.

## Repositorio compartido y entrega a mis compañeros

Organiza el proyecto para que la lógica reutilizable esté separada de los notebooks de resultados. Como mínimo, incluye el simulador corregido, las funciones teóricas, el modelo SIR, mi notebook, el anexo del código original, un `README` con instrucciones de ejecución y las dependencias. Documenta las entradas, salidas, unidades y supuestos de las funciones.

Gerardo necesitará reutilizar el simulador al cambiar distribuciones de estancia y capacidad de camas. Juan Diego necesitará probar distintos números de camas y, posteriormente, llegadas variables basadas en la curva SIR. Deja claras las **extensiones necesarias** sin presentar como implementadas funciones que aún no existan. Mantén el código original de IA intacto y separado del módulo que usarán los demás.

Al terminar, dame: **(a)** una lista de requisitos de mis incisos cumplidos con el archivo donde aparece cada uno; **(b)** los archivos concretos que debo pasar a Gerardo y Juan Diego; y **(c)** cualquier verificación que yo deba hacer personalmente antes de entregar y defender el trabajo. Prioriza resultados verificables y explicaciones que pueda comprender, corregir y defender.

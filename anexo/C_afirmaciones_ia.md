# Anexo C — Afirmaciones textuales de la IA sobre la corrección del código original (Task 1.2b)

Estas son las únicas afirmaciones que la IA (Claude Sonnet 5.5, en Claude Code) hizo sobre el código de
`B_simulador_original_ia.py` **en el momento de generarlo y archivarlo**, antes de ejecutarlo o revisarlo.
Se copian literalmente.

## C.1 Afirmaciones dentro del propio archivo (docstring y comentarios)

1. `"Primera version generada por IA (Claude Sonnet 5.5) en Claude Code."`
2. `"Usa una lista de eventos futuros (heap). Tiempos en dias."`
3. Comentario: `"# cerrar el ultimo tramo hasta el horizonte"`
4. Comentario: `"# acumular areas del tramo anterior (solo despues del calentamiento)"`
5. Docstring de `replicate`: `"Corre n_reps replicas y devuelve media e IC 95% por medida."`

## C.2 Afirmaciones en la conversación al archivar el código

Mensaje literal de la IA inmediatamente antes de escribir el archivo:
> "Ya tengo los requisitos de la hoja. Sigo con tu orden de trabajo: archivo primero el prompt y la primera
> versión del simulador, sin ejecutarla ni revisarla."

Mensaje literal tras escribirlo: ninguna afirmación de corrección. La IA **no** afirmó que el código funcione,
coincida con la teoría ni haya pasado pruebas, y de hecho **no lo había ejecutado** al archivarlo.

## C.3 Implicación para el Task 1.4c

Como la IA no hizo afirmaciones de corrección más allá de las del punto C.1 (afirmaciones *implícitas* del código:
p. ej. "cierra el último tramo hasta el horizonte", "IC 95 % por medida"), el ejercicio 1.4c se aplica a esas
afirmaciones implícitas (ver el notebook `notebooks/task1_teoria_simulador_y_validacion.ipynb` y el informe del grupo).
Si la IA hiciera afirmaciones adicionales en la conversación sobre este código, **el autor debe copiarlas aquí
literalmente** (no se deben reconstruir de memoria).

## C.4 Integridad

`SHA256SUMS.txt` contiene el hash de `B_simulador_original_ia.py` calculado inmediatamente después de escribirlo y
antes de ejecutarlo. Se comprueba con `cd anexo && shasum -a 256 -c SHA256SUMS.txt` y también en la sección Task 1.2
del notebook de Task 1. Nota: el hash de `A_prompt_exacto.md` corresponde a la transcripción, no necesariamente al
texto original enviado.

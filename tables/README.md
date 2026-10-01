# Tablas del TFM — Procedencia y notas

Este directorio contiene las tablas del trabajo. Las tablas se dividen en dos
grupos conorigine muy distinta.

## 1. Tablas reproducibles (generadas por `R/analysis.R`)

Estas se recalculan desde cero a partir de `data/analysis_data.csv`:

- `clinical_descriptives.csv`
- `pca_loadings.csv`
- `cluster_clinical_comparison.csv`

Para regenerarlas:

```bash
Rscript R/analysis.R
```

## 2. Agregados históricos estáticos (NO reproducibles)

Las siguientes dos tablas son una transcripción literal de la tesis ya
entregada. **No** se reconstruyen a partir de datos de nivel participante:

- `sample_characteristics.csv` → Tabla 1 de la tesis (PDF, págs. 26-27)
  «Características sociodemográficas y contextuales de la muestra»
- `cluster_sociodemographic_comparison.csv` → Tabla 4 de la tesis (PDF, pág. 40)
  «Comparación sociodemográfica entre los perfiles clínicos»

Fuente: `TFM Pablo Pérez Vargas.pdf`.

### Por qué son estáticas

El conjunto de datos público (`data/analysis_data.csv`) excluye deliberadamente
las variables sociodemográficas (edad exacta, género, nivel de estudios,
situación laboral, tamaño de empresa, jornada, situación económica, residencia,
composición del hogar, apoyo percibido). Por tanto, **estas tablas no pueden
volver a calcularse** a partir de los datos publicados y se conservan aquí como
registro histórico de los resultados tal como se entregaron.

Tampoco se incluyen análisis de fiabilidad (alfa de Cronbach), ya que el
conjunto público no incluye los ítems de cuestionario.

### Nota sobre la edad

En la Tabla 1 del texto original la edad se describe en el cuerpo del texto
(`M = 28.5`, `DT = 8.88`) y **no** aparece como fila de la tabla; por eso no se
incluye en `sample_characteristics.csv`. Sí aparece en la Tabla 4 como
`Edad (Media, DT)`, que sí se ha transcrito.

### Formato

Ambos CSV comienzan por una fila de encabezado real (`Variable,Category,...`),
sin líneas de comentario, para que se puedan leer con cualquier lector CSV
estándar. Esta nota es la que documenta la procedencia.

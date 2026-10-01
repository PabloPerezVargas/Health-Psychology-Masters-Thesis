# ==============================================================================
# create_analysis_dataset.R
#
# Builds a de-identified analysis dataset from the raw thesis data.
#
#   Input  : R TFM/datos91.csv   (read-only; never modified)
#   Output : data/analysis_data.csv
#
# Scoring logic is reproduced from "R TFM/R SCRIPT TFM MPGS.R" (section 3),
# including the reverse-scored items and the IDERE formula.
#
# The output contains ONLY an anonymous sequential id, one contextual count
# and the ten derived instrument scores. No raw items, no demographics,
# no timestamps and no free text are exported.
#
# Run with:  Rscript create_analysis_dataset.R
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

# ------------------------------------------------------------------------------
# 0. PATHS AND SAFETY GUARDS
# ------------------------------------------------------------------------------

obtener_dir_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  coincidencia <- grep("^--file=", args, value = TRUE)
  if (length(coincidencia) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", coincidencia[1]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

DIR_BASE     <- obtener_dir_script()
RUTA_ENTRADA <- file.path(DIR_BASE, "R TFM", "datos91.csv")
DIR_SALIDA   <- file.path(DIR_BASE, "data")
RUTA_SALIDA  <- file.path(DIR_SALIDA, "analysis_data.csv")

if (!file.exists(RUTA_ENTRADA)) {
  stop("Input file not found: ", RUTA_ENTRADA)
}

# Never overwrite an existing output file.
if (file.exists(RUTA_SALIDA)) {
  stop("Output already exists and will NOT be overwritten: ", RUTA_SALIDA)
}

if (!dir.exists(DIR_SALIDA)) {
  dir.create(DIR_SALIDA, recursive = TRUE)
  cat("Created folder:", DIR_SALIDA, "\n")
}

# ------------------------------------------------------------------------------
# 1. READ AND PREPARE  (same preparation as the original script)
# ------------------------------------------------------------------------------

datos <- read.csv(RUTA_ENTRADA, fileEncoding = "UTF-8") %>%
  dplyr::select(-(1:2)) %>%
  dplyr::mutate(id = row_number()) %>%
  dplyr::relocate(id)

n_filas_crudo <- nrow(datos)

names(datos)[1:14] <- c(
  "id", "Edad", "Genero", "NivelEstudios", "SituacionLaboral",
  "SectorLaboral", "TamanoEmpresa", "JornadaLaboral", "SituacionEconomica",
  "Residencia", "ConQuienVive", "ResponsabilidadCuidado",
  "EventosEstresantes", "NivelApoyo"
)

colnames(datos)[15:42] <- paste0("CTQ_", 1:28)
colnames(datos)[43:82] <- paste0("STAI_", 1:40)
colnames(datos)[83:124] <- paste0("IDERE_", 1:42)

# ------------------------------------------------------------------------------
# 2. SAMPLE FILTER AND CONTEXTUAL VARIABLE
# ------------------------------------------------------------------------------

datos <- datos %>%
  mutate(Edad = as.numeric(gsub("[^0-9]", "", Edad))) %>%
  filter(Edad >= 18 & Edad <= 70) %>%
  mutate(
    EventosEstresantes = sapply(EventosEstresantes, function(x) {
      if (is.na(x) || x == "" || grepl("Ninguno de los anteriores", x)) return(0)
      length(unlist(strsplit(x, ",")))
    })
  )

n_filas_finales <- nrow(datos)
cat("Rows read:", n_filas_crudo, "| Rows after age filter (18-70):", n_filas_finales, "\n")

# ------------------------------------------------------------------------------
# 3. SCORING  (identical to section 3 of "R SCRIPT TFM MPGS.R")
# ------------------------------------------------------------------------------

# A) CTQ - reverse-score the designated items first, then sum subscales.
ctq_itemsinversos <- paste0("CTQ_", c(2, 5, 7, 13, 19, 26, 28))
datos <- datos %>%
  mutate(across(all_of(ctq_itemsinversos), ~ 6 - .)) %>%
  mutate(
    CTQ_AE = rowSums(pick(paste0("CTQ_", c(3, 8, 14, 18, 25))), na.rm = TRUE),
    CTQ_AF = rowSums(pick(paste0("CTQ_", c(9, 11, 12, 15, 17))), na.rm = TRUE),
    CTQ_AS = rowSums(pick(paste0("CTQ_", c(20, 21, 23, 24, 27))), na.rm = TRUE),
    CTQ_NE = rowSums(pick(paste0("CTQ_", c(5, 7, 13, 19, 28))), na.rm = TRUE),
    CTQ_NF = rowSums(pick(paste0("CTQ_", c(1, 2, 4, 6, 26))), na.rm = TRUE),
    CTQ_TOTAL = CTQ_AE + CTQ_AF + CTQ_AS + CTQ_NE + CTQ_NF
  )

# B) STAI - reverse-score the designated items, then sum state and trait.
stai_itemsinversos <- paste0("STAI_", c(1, 2, 5, 8, 10, 11, 15, 16, 19,
                                        20, 21, 26, 27, 30, 33, 36, 39))
datos <- datos %>%
  mutate(across(all_of(stai_itemsinversos), ~ 5 - .)) %>%
  mutate(
    STAI_State = rowSums(pick(paste0("STAI_", 1:20)), na.rm = TRUE),
    STAI_Trait = rowSums(pick(paste0("STAI_", 21:40)), na.rm = TRUE)
  )

# C) IDERE - direct formula, no column inversion.
idere_est_A <- paste0("IDERE_", c(1, 2, 6, 8, 9, 10, 13, 15, 18, 20))
idere_est_B <- paste0("IDERE_", c(3, 4, 5, 7, 11, 12, 14, 16, 17, 19))
idere_ras_A <- paste0("IDERE_", c(21, 23, 25, 27, 28, 31, 34, 36, 37, 40, 42))
idere_ras_B <- paste0("IDERE_", c(22, 24, 26, 29, 30, 32, 33, 35, 38, 39, 41))

datos <- datos %>%
  mutate(
    IDERE_State = rowSums(pick(all_of(idere_est_A)), na.rm = TRUE) -
      rowSums(pick(all_of(idere_est_B)), na.rm = TRUE) + 50,
    IDERE_Trait = rowSums(pick(all_of(idere_ras_A)), na.rm = TRUE) -
      rowSums(pick(all_of(idere_ras_B)), na.rm = TRUE) + 50
  )

# ------------------------------------------------------------------------------
# 4. BUILD THE OUTPUT DATA FRAME (12 columns, nothing else)
# ------------------------------------------------------------------------------

analysis_data <- data.frame(
  participant_id     = seq_len(n_filas_finales),
  stressful_events_n = as.integer(datos$EventosEstresantes),
  CTQ_AE    = as.numeric(datos$CTQ_AE),
  CTQ_AF    = as.numeric(datos$CTQ_AF),
  CTQ_AS    = as.numeric(datos$CTQ_AS),
  CTQ_NE    = as.numeric(datos$CTQ_NE),
  CTQ_NF    = as.numeric(datos$CTQ_NF),
  CTQ_TOTAL = as.numeric(datos$CTQ_TOTAL),
  STAI_State = as.numeric(datos$STAI_State),
  STAI_Trait = as.numeric(datos$STAI_Trait),
  IDERE_State = as.numeric(datos$IDERE_State),
  IDERE_Trait = as.numeric(datos$IDERE_Trait),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 5. VERIFICATION (runs before the file is written)
# ------------------------------------------------------------------------------

columnas_esperadas <- c(
  "participant_id", "stressful_events_n",
  "CTQ_AE", "CTQ_AF", "CTQ_AS", "CTQ_NE", "CTQ_NF", "CTQ_TOTAL",
  "STAI_State", "STAI_Trait", "IDERE_State", "IDERE_Trait"
)

errores <- character(0)

# 5.1 Dimensions and column names -------------------------------------------
if (nrow(analysis_data) != n_filas_finales) errores <- c(errores, "unexpected row count")
if (!identical(names(analysis_data), columnas_esperadas)) errores <- c(errores, "column names do not match")

# 5.2 No raw item-level or identifying variables ---------------------------
patron_item <- "^(CTQ|STAI|IDERE)_[0-9]+$"
items_encontrados <- grep(patron_item, names(analysis_data), value = TRUE)
if (length(items_encontrados) > 0) errores <- c(errores, "raw item-level columns present")

prohibidas <- c("Edad", "Genero", "NivelEstudios", "SituacionLaboral", "SectorLaboral",
                "TamanoEmpresa", "JornadaLaboral", "SituacionEconomica", "Residencia",
                "ConQuienVive", "ResponsabilidadCuidado", "NivelApoyo", "id",
                "Marca.temporal")
presentes <- intersect(prohibidas, names(analysis_data))
if (length(presentes) > 0) errores <- c(errores, "prohibited raw columns present")

# 5.3 Missing values -------------------------------------------------------
n_missing <- sum(is.na(analysis_data))
if (n_missing > 0) errores <- c(errores, "missing values found")

# 5.4 Independent recomputation of the scores from the raw CSV -------------
ref <- read.csv(RUTA_ENTRADA, fileEncoding = "UTF-8")
ref$edad_num <- as.numeric(gsub("[^0-9]", "", ref[[3]]))
ref <- ref[!is.na(ref$edad_num) & ref$edad_num >= 18 & ref$edad_num <= 70, ]

ctq_ref <- as.matrix(ref[, 16:43])
for (j in c(2, 5, 7, 13, 19, 26, 28)) ctq_ref[, j] <- 6 - ctq_ref[, j]
ref$CTQ_AE <- rowSums(ctq_ref[, c(3, 8, 14, 18, 25)], na.rm = TRUE)
ref$CTQ_AF <- rowSums(ctq_ref[, c(9, 11, 12, 15, 17)], na.rm = TRUE)
ref$CTQ_AS <- rowSums(ctq_ref[, c(20, 21, 23, 24, 27)], na.rm = TRUE)
ref$CTQ_NE <- rowSums(ctq_ref[, c(5, 7, 13, 19, 28)], na.rm = TRUE)
ref$CTQ_NF <- rowSums(ctq_ref[, c(1, 2, 4, 6, 26)], na.rm = TRUE)
ref$CTQ_TOTAL <- ref$CTQ_AE + ref$CTQ_AF + ref$CTQ_AS + ref$CTQ_NE + ref$CTQ_NF

stai_ref <- as.matrix(ref[, 44:83])
for (j in c(1, 2, 5, 8, 10, 11, 15, 16, 19, 20, 21, 26, 27, 30, 33, 36, 39)) stai_ref[, j] <- 5 - stai_ref[, j]
ref$STAI_State <- rowSums(stai_ref[, 1:20], na.rm = TRUE)
ref$STAI_Trait <- rowSums(stai_ref[, 21:40], na.rm = TRUE)

idere_ref <- as.matrix(ref[, 84:125])
ref$IDERE_State <- rowSums(idere_ref[, c(1, 2, 6, 8, 9, 10, 13, 15, 18, 20)], na.rm = TRUE) -
                   rowSums(idere_ref[, c(3, 4, 5, 7, 11, 12, 14, 16, 17, 19)], na.rm = TRUE) + 50
ref$IDERE_Trait <- rowSums(idere_ref[, c(21, 23, 25, 27, 28, 31, 34, 36, 37, 40, 42)], na.rm = TRUE) -
                   rowSums(idere_ref[, c(22, 24, 26, 29, 30, 32, 33, 35, 38, 39, 41)], na.rm = TRUE) + 50

puntuaciones <- c("CTQ_AE", "CTQ_AF", "CTQ_AS", "CTQ_NE", "CTQ_NF", "CTQ_TOTAL",
                  "STAI_State", "STAI_Trait", "IDERE_State", "IDERE_Trait")

coincidencias <- vapply(puntuaciones, function(v) {
  isTRUE(all.equal(as.numeric(analysis_data[[v]]), as.numeric(ref[[v]]), tolerance = 1e-9))
}, logical(1))
if (any(!coincidencias)) {
  errores <- c(errores, paste("scores differ from reference:", paste(names(coincidencias)[!coincidencias], collapse = ", ")))
}

# 5.5 Instrument ranges ---------------------------------------------------
rangos <- list(CTQ_AE = c(5, 25), CTQ_AF = c(5, 25), CTQ_AS = c(5, 25),
               CTQ_NE = c(5, 25), CTQ_NF = c(5, 25), CTQ_TOTAL = c(25, 125),
               STAI_State = c(20, 80), STAI_Trait = c(20, 80),
               IDERE_State = c(20, 80), IDERE_Trait = c(20, 80))
fuera <- 0
for (v in names(rangos)) {
  fuera <- fuera + sum(analysis_data[[v]] < rangos[[v]][1] |
                       analysis_data[[v]] > rangos[[v]][2], na.rm = TRUE)
}
if (fuera > 0) errores <- c(errores, "values outside instrument range")

# ------------------------------------------------------------------------------
# 6. REPORT
# ------------------------------------------------------------------------------

cat("\n================ VERIFICATION ================\n")
cat("Rows:", nrow(analysis_data), "\n")
cat("Columns:", ncol(analysis_data), "\n")
cat("Column names:\n")
cat(paste0("  ", seq_along(names(analysis_data)), ". ", names(analysis_data)), sep = "\n")
cat("\nRaw item-level columns present:",
    if (length(items_encontrados) == 0) "NONE" else paste(items_encontrados, collapse = ", "), "\n")
cat("Prohibited raw/demographic columns present:",
    if (length(presentes) == 0) "NONE" else paste(presentes, collapse = ", "), "\n")
cat("Total NA values:", n_missing, "\n")
cat("Values outside instrument range:", fuera, "\n")
cat("\nScore verification (independent recomputation):\n")
cat(paste0("  ", names(coincidencias), ": ",
           ifelse(coincidencias, "MATCH", "MISMATCH")), sep = "\n")
cat("\nDescriptive summary:\n")
print(summary(analysis_data))

# ------------------------------------------------------------------------------
# 7. WRITE (only if every check passed)
# ------------------------------------------------------------------------------

if (length(errores) > 0) {
  stop("\nVERIFICATION FAILED. No file was written:\n  - ",
       paste(errores, collapse = "\n  - "))
}

write.csv(analysis_data, RUTA_SALIDA, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n============================================\n")
cat("Dataset written to:", RUTA_SALIDA, "\n")
cat("Rows:", nrow(analysis_data), "| Columns:", ncol(analysis_data), "\n")
cat("============================================\n")

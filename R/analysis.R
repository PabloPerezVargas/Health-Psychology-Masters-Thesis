# ==============================================================================
# R/analysis.R
#
# Main analyses for the thesis, reproduced from "R TFM/R SCRIPT TFM MPGS.R".
#
#   Input  : data/analysis_data.csv  (only dataset used)
#   Outputs: figures/pca_biplot.png
#            figures/cluster_kmeans.png
#            figures/silhouette.png
#            console report (descriptives, KMO, Bartlett, PCA, clustering)
#
# Paths are resolved relative to the repository root, derived from the location
# of this script, so the analysis runs from anywhere without setwd().
#
# Run with:  Rscript R/analysis.R
# ==============================================================================

suppressPackageStartupMessages({
  library(psych)
  library(cluster)
  library(factoextra)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------------------------

# Repository root = parent directory of the folder containing this script.
obtener_dir_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  coincidencia <- grep("^--file=", args, value = TRUE)
  if (length(coincidencia) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", coincidencia[1]), winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

DIR_REPO     <- normalizePath(file.path(obtener_dir_script(), ".."), winslash = "/")
RUTA_DATOS   <- file.path(DIR_REPO, "data", "analysis_data.csv")
DIR_FIGURAS  <- file.path(DIR_REPO, "figures")
DIR_TABLES   <- file.path(DIR_REPO, "tables")

if (!file.exists(RUTA_DATOS)) {
  stop("Input dataset not found: ", RUTA_DATOS)
}

# Create figures/ automatically if it does not exist.
if (!dir.exists(DIR_FIGURAS)) {
  dir.create(DIR_FIGURAS, recursive = TRUE)
  cat("Created figures folder:", DIR_FIGURAS, "\n")
}

# Create tables/ automatically if it does not exist.
if (!dir.exists(DIR_TABLES)) {
  dir.create(DIR_TABLES, recursive = TRUE)
  cat("Created tables folder:", DIR_TABLES, "\n")
}

# ------------------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------------------

datos <- read.csv(RUTA_DATOS, fileEncoding = "UTF-8")

# Analysis variables, matching the composition used in the original script:
# number of stressful events + 5 CTQ subscales + STAI + IDERE.
vars_analisis <- c(
  "stressful_events_n",
  "CTQ_AE", "CTQ_AF", "CTQ_AS", "CTQ_NE", "CTQ_NF",
  "STAI_State", "STAI_Trait", "IDERE_State", "IDERE_Trait"
)

cat("Rows:", nrow(datos), "| Columns:", ncol(datos), "\n\n")

# ------------------------------------------------------------------------------
# 2. DESCRIPTIVE STATISTICS FOR THE CLINICAL VARIABLES
#    (n, mean, standard deviation, min, max - as in section 4.2 of the original)
# ------------------------------------------------------------------------------

calc_stats <- function(df, nivel_sangria = 0) {
  data.frame(
    Variable   = names(df),
    n          = sapply(df, function(x) sum(!is.na(x))),
    M          = sapply(df, function(x) round(mean(x, na.rm = TRUE), 2)),
    DT         = sapply(df, function(x) round(sd(x, na.rm = TRUE), 2)),
    Min        = sapply(df, function(x) min(x, na.rm = TRUE)),
    Max        = sapply(df, function(x) max(x, na.rm = TRUE)),
    indent     = nivel_sangria,
    is_header  = FALSE,
    stringsAsFactors = FALSE
  )
}

vars_simples   <- datos %>% select("Eventos Estresantes (Nº)" = stressful_events_n)
vars_trauma    <- datos %>% select("Trauma Infantil (CTQ Total)" = CTQ_TOTAL,
                                    "Abuso Emocional" = CTQ_AE,
                                    "Abuso Físico" = CTQ_AF,
                                    "Abuso Sexual" = CTQ_AS,
                                    "Negligencia Emocional" = CTQ_NE,
                                    "Negligencia Física" = CTQ_NF)
vars_ansiedad  <- datos %>% select("Estado" = STAI_State, "Rasgo" = STAI_Trait)
vars_depresion <- datos %>% select("Estado" = IDERE_State, "Rasgo" = IDERE_Trait)

tabla_descriptivos <- rbind(
  calc_stats(vars_simples, 0),
  calc_stats(vars_trauma[1], 0), calc_stats(vars_trauma[-1], 1),
  data.frame(Variable = "Ansiedad (STAI)", n = NA, M = NA, DT = NA, Min = NA, Max = NA,
             indent = 0, is_header = TRUE, stringsAsFactors = FALSE),
  calc_stats(vars_ansiedad, 1),
  data.frame(Variable = "Depresión (IDERE)", n = NA, M = NA, DT = NA, Min = NA, Max = NA,
             indent = 0, is_header = TRUE, stringsAsFactors = FALSE),
  calc_stats(vars_depresion, 1)
)

cat("=== 1. ESTADÍSTICOS DESCRIPTIVOS (VARIABLES CLÍNICAS) ===\n")
print(tabla_descriptivos[, c("Variable", "n", "M", "DT", "Min", "Max")], row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. FIGURE 1 - DISTRIBUTION OF THE FIVE CTQ SUBSCALES (FLOOR EFFECT)
#     Reproduces the "efecto suelo" boxplot of the original thesis: one boxplot
#     per subscale, with the mean overlaid as a white diamond.
# ------------------------------------------------------------------------------

cat("\n=== 2. DISTRIBUCIÓN DE LAS SUBESCALAS DEL CTQ ===\n")

datos_long_ctq <- datos %>%
  select(CTQ_AE, CTQ_AF, CTQ_AS, CTQ_NE, CTQ_NF) %>%
  pivot_longer(cols = everything(),
               names_to = "Subescala", values_to = "Puntuacion")

# Spanish labels matching the terminology of the thesis.
etiquetas_ctq <- c(
  CTQ_AE = "Abuso Emocional",
  CTQ_AF = "Abuso Físico",
  CTQ_AS = "Abuso Sexual",
  CTQ_NE = "Negligencia Emocional",
  CTQ_NF = "Negligencia Física"
)
datos_long_ctq$Subescala <- factor(datos_long_ctq$Subescala,
                                   levels = names(etiquetas_ctq),
                                   labels = unname(etiquetas_ctq))

figura_ctq <- ggplot(datos_long_ctq, aes(x = Subescala, y = Puntuacion, fill = Subescala)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white", size = 2) +
  labs(title = "Distribución de Trauma (Efecto Suelo)",
       x = NULL, y = "Puntuación CTQ") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
print(figura_ctq)

ruta_ctq <- file.path(DIR_FIGURAS, "ctq_subscale_distributions.png")
ggsave(filename = ruta_ctq, plot = figura_ctq, width = 9, height = 6, dpi = 300)
cat("CTQ subscale distributions saved to:", ruta_ctq, "\n")

# ------------------------------------------------------------------------------
# 3. KMO AND BARTLETT'S TEST OF SPHERICITY
#    (applied to the analysis matrix, excluding the identifier column)
# ------------------------------------------------------------------------------

matriz_analisis <- datos %>% select(all_of(vars_analisis))

cat("\n=== 3. KMO Y TEST DE BARTLETT ===\n")
kmo_resultado <- psych::KMO(matriz_analisis)
cat("KMO (overall):", round(kmo_resultado$MSA, 3), "\n\n")
print(round(kmo_resultado$MSAi, 3))

bartlett_resultado <- psych::cortest.bartlett(matriz_analisis)
cat("\nBartlett's test of sphericity\n")
cat("  Chi-squared =", round(bartlett_resultado$chisq, 2), "\n")
cat("  df          =", bartlett_resultado$df, "\n")
cat("  p-value     =", format.pval(bartlett_resultado$p.value, eps = 0.001), "\n")

# ------------------------------------------------------------------------------
# 4. PCA (prcomp, standardised) AND PC1 SIGN ORIENTATION
#    The original flips the sign of PC1 in both the rotation and the scores.
# ------------------------------------------------------------------------------

cat("\n=== 4. ANÁLISIS DE COMPONENTES PRINCIPALES (PCA) ===\n")
pca_res <- prcomp(as.matrix(matriz_analisis), scale. = TRUE)

# Same PC1 sign orientation as the original script.
pca_res$rotation[, 1] <- -1 * pca_res$rotation[, 1]
pca_res$x[, 1]        <- -1 * pca_res$x[, 1]

varianza_explicada <- pca_res$sdev^2 / sum(pca_res$sdev^2)
etiquetas <- sprintf("PC%d", seq_along(varianza_explicada))

print(round(varianza_explicada * 100, 2))
cat("\n% of variance explained by PC1:",
    round(varianza_explicada[1] * 100, 2), "%\n")
cat("% of variance explained by PC2:",
    round(varianza_explicada[2] * 100, 2), "%\n")
cat("Cumulative variance (PC1 + PC2):",
    round(sum(varianza_explicada[1:2]) * 100, 2), "%\n")

# Contributions of each variable to PC1 and PC2.
contrib <- fviz_contrib(pca_res, choice = "var", axes = 1, top = 10, fill = "steelblue")
print(contrib)
ggsave(filename = file.path(DIR_FIGURAS, "pca_contrib_pc1.png"),
       plot = contrib, width = 8, height = 6, dpi = 300)

contrib2 <- fviz_contrib(pca_res, choice = "var", axes = 2, top = 10, fill = "darkred")
print(contrib2)
ggsave(filename = file.path(DIR_FIGURAS, "pca_contrib_pc2.png"),
       plot = contrib2, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5. PCA LOADINGS
# ------------------------------------------------------------------------------

cat("\n=== 5. PESAS (LOADINGS) DE LOS COMPONENTES ===\n")
print(round(pca_res$rotation, 3))

cat("\nLoadings on PC1 and PC2 (variables with |loading| >= 0.30):\n")
loadings_2d <- data.frame(
  Variable    = rownames(pca_res$rotation),
  PC1         = pca_res$rotation[, 1],
  PC2         = pca_res$rotation[, 2],
  row.names   = NULL
)
print(loadings_2d[abs(loadings_2d$PC1) >= 0.30 | abs(loadings_2d$PC2) >= 0.30, ])

# ------------------------------------------------------------------------------
# 6. PCA BIPLOT -> figures/pca_biplot.png
# ------------------------------------------------------------------------------

biplot <- fviz_pca_biplot(pca_res,
                          label = "var", repel = TRUE,
                          col.var = "black", col.ind = "grey40",
                          alpha.ind = 0.6, pointsize = 1.5,
                          title = "Biplot PCA")
print(biplot)

ruta_biplot <- file.path(DIR_FIGURAS, "pca_biplot.png")
ggsave(filename = ruta_biplot, plot = biplot, width = 12, height = 7, dpi = 300)
cat("\nBiplot saved to:", ruta_biplot, "\n")

# ------------------------------------------------------------------------------
# 7. CLUSTERING TENDENCY (HOPKINS) AND SILHOUETTE METHOD
# ------------------------------------------------------------------------------

# Scores of the first three components, as in the original script.
pca_scores <- as.data.frame(pca_res$x[, 1:3])

cat("\n=== 6. TENDENCIA AL AGRUPAMIENTO (HOPKINS) ===\n")
res_hopkins <- get_clust_tendency(pca_scores, n = nrow(pca_scores) - 1, graph = FALSE)
cat("Hopkins statistic:", round(res_hopkins$hopkins_stat, 4), "\n")

cat("\n=== 7. MÉTODO DE LA SILUETA ===\n")
set.seed(123)
grafico_silueta <- fviz_nbclust(pca_scores, kmeans, method = "silhouette") +
  labs(title = "Método de la Silueta")
print(grafico_silueta)
ggsave(filename = file.path(DIR_FIGURAS, "silhouette.png"),
       plot = grafico_silueta, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 8. K-MEANS CLUSTERING (K = 2, seed 123)
# ------------------------------------------------------------------------------

cat("\n=== 8. CLUSTERING K-MEANS (K = 2) ===\n")
set.seed(123)
kmeans_final <- kmeans(pca_scores, centers = 2, nstart = 25)

cat("Cluster sizes:\n")
print(kmeans_final$size)

datos$Cluster <- as.factor(kmeans_final$cluster)

cat("\nTotal within-cluster sum of squares (tot.withinss):",
    round(kmeans_final$tot.withinss, 2), "\n")

grafico_cluster <- fviz_cluster(kmeans_final,
                               data = pca_scores[, 1:2],
                               palette = c("#2E9FDF", "#E7B800"),
                               geom = "point", ellipse.type = "norm",
                               ellipse.level = 0.95, ggtheme = theme_minimal(),
                               main = "Clustering K-Means (K=2)")
print(grafico_cluster)
ggsave(filename = file.path(DIR_FIGURAS, "cluster_kmeans.png"),
       plot = grafico_cluster, width = 9, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 9. CLUSTER CLINICAL PROFILES
# ------------------------------------------------------------------------------

cat("\n=== 9. PERFIL CLÍNICO POR CLÚSTER ===\n")
tabla_clinica <- datos %>%
  group_by(Cluster) %>%
  summarise(N        = n(),
            Trauma   = mean(CTQ_TOTAL),
            Ansiedad = mean(STAI_Trait),
            Depresion = mean(IDERE_Trait),
            Eventos  = mean(stressful_events_n),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), round, 2))

print(tabla_clinica)

# ------------------------------------------------------------------------------
# 10. BETWEEN-CLUSTER COMPARISONS (independent-samples t-tests)
#     Same four variables as the original script.
# ------------------------------------------------------------------------------

cat("\n=== 10. DIFERENCIAS CLÍNICAS ENTRE CLÚSTERES (T-STUDENT) ===\n")
vars_clinicas <- c("CTQ_TOTAL", "STAI_Trait", "IDERE_Trait", "stressful_events_n")

for (variable in vars_clinicas) {
  test <- t.test(datos[[variable]] ~ datos$Cluster)
  significacion <- ifelse(test$p.value < 0.05, "Significant", "Not significant")
  cat(sprintf("%-20s | t = %7.3f | p = %s -> %s\n",
              variable,
              unname(test$statistic),
              format.pval(test$p.value, eps = 0.001),
              significacion))
}

# ------------------------------------------------------------------------------
# 11. REPRODUCIBLE TABLES (CSV) -> tables/
#     Only clinical/contextual variables: the public dataset contains no
#     sociodemographic data and no item-level responses, so neither participant
#     profile tables nor reliability (Cronbach's alpha) tables are produced.
# ------------------------------------------------------------------------------

cat("\n=== 11. TABLAS REPRODUCIBLES (CSV) ===\n")

# 11.1 Clinical descriptives ------------------------------------------------
vars_descriptivas <- c("stressful_events_n", "CTQ_TOTAL",
                       "CTQ_AE", "CTQ_AF", "CTQ_AS", "CTQ_NE", "CTQ_NF",
                       "STAI_State", "STAI_Trait", "IDERE_State", "IDERE_Trait")

tabla_clinica_descriptivos <- data.frame(
  variable = vars_descriptivas,
  n        = sapply(datos[vars_descriptivas], function(x) sum(!is.na(x))),
  mean     = round(sapply(datos[vars_descriptivas], function(x) mean(x, na.rm = TRUE)), 2),
  sd       = round(sapply(datos[vars_descriptivas], function(x) sd(x, na.rm = TRUE)), 2),
  min      = sapply(datos[vars_descriptivas], function(x) min(x, na.rm = TRUE)),
  max      = sapply(datos[vars_descriptivas], function(x) max(x, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

ruta_tabla_descriptivos <- file.path(DIR_TABLES, "clinical_descriptives.csv")
write.csv(tabla_clinica_descriptivos, ruta_tabla_descriptivos,
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Written:", ruta_tabla_descriptivos, "\n")
print(tabla_clinica_descriptivos, row.names = FALSE)

# 11.2 PCA loadings (PC1 and PC2, all variables used in the PCA) -----------
tabla_loadings <- data.frame(
  variable = rownames(pca_res$rotation),
  PC1      = round(pca_res$rotation[, 1], 6),
  PC2      = round(pca_res$rotation[, 2], 6),
  pct_variance_PC1 = round(varianza_explicada[1] * 100, 2),
  pct_variance_PC2 = round(varianza_explicada[2] * 100, 2),
  stringsAsFactors = FALSE
)

ruta_tabla_loadings <- file.path(DIR_TABLES, "pca_loadings.csv")
write.csv(tabla_loadings, ruta_tabla_loadings,
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nWritten:", ruta_tabla_loadings, "\n")
print(tabla_loadings, row.names = FALSE)

# 11.3 Clinical cluster comparison (n, mean, SD and the t-test p-value) -----
# One row per variable and cluster; the p-value comes from the same
# independent-samples t-test used in the console comparison above.
filas_cluster <- list()
for (variable in vars_clinicas) {
  test <- t.test(datos[[variable]] ~ datos$Cluster)

  for (cluster in levels(datos$Cluster)) {
    valores <- datos[[variable]][datos$Cluster == cluster]
    filas_cluster[[length(filas_cluster) + 1]] <- data.frame(
      variable = variable,
      cluster  = cluster,
      n        = length(valores),
      mean     = round(mean(valores, na.rm = TRUE), 2),
      sd       = round(sd(valores, na.rm = TRUE), 2),
      p_value  = signif(test$p.value, 6),
      stringsAsFactors = FALSE
    )
  }
}
tabla_cluster_comparacion <- do.call(rbind, filas_cluster)

ruta_tabla_cluster <- file.path(DIR_TABLES, "cluster_clinical_comparison.csv")
write.csv(tabla_cluster_comparacion, ruta_tabla_cluster,
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nWritten:", ruta_tabla_cluster, "\n")
print(tabla_cluster_comparacion, row.names = FALSE)

# ------------------------------------------------------------------------------
# 12. DONE
# ------------------------------------------------------------------------------

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Figures written to:", DIR_FIGURAS, "\n")
cat("  - ctq_subscale_distributions.png\n")
cat("  - pca_contrib_pc1.png\n")
cat("  - pca_contrib_pc2.png\n")
cat("  - pca_biplot.png\n")
cat("  - silhouette.png\n")
cat("  - cluster_kmeans.png\n")
cat("Tables written to:", DIR_TABLES, "\n")
cat("  - clinical_descriptives.csv\n")
cat("  - pca_loadings.csv\n")
cat("  - cluster_clinical_comparison.csv\n")

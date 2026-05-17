################################################################################
#
#   title: "Gene-level differential expression analysis using DESeq2"
#   author: "Verónica Cabeza de Vaca Tocino"
#   date: "2026-04"
#
################################################################################

## 1. Setup y Carga de librerías ####

library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)

getwd()

# Crear directorio para guardar resultados si no existe
dataset_name <- "GSE164541"

output_dir <- file.path("Output", dataset_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Definir un tema estándar para todas las gráficas (Estilo Paper)
tema_paper <- theme_classic(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )

## 2. Carga y Preparación de Datos ####

# Carga de la matriz de conteos
data_raw <- read.csv("Data/GSE164541/GSE164541_ANT_count.csv", header = TRUE, row.names = 1)
head(data_raw)

# Ordenar/filtrar columnas para que cuadre con la metadata
data <- data_raw[3:17] %>%
  dplyr::select(PN1, PN2, PN3, PN4, PN5, PA1, PA2, PA3, PA4, PA5, PT1, PT2, PT3, PT4, PT5)
head(data)
dim(data)

# Carga de metadata
meta <- read.csv("Data/GSE164541/GSE164541_filtered_metadata.csv", header = TRUE) %>%
  data.frame(row.names = colnames(data))
head(meta)

# Generar metadata (coldata) extrayendo información de los nombres de muestra
coldata <- data.frame(
  Patient = as.factor(substr(colnames(data), 3, 3)),
  Condition = as.factor(substr(colnames(data), 1, 2)),
  row.names = colnames(data)
)

# Relevel para que Normal (PN) sea la referencia basal
coldata$Condition <- relevel(coldata$Condition, ref = "PN")

# Check de seguridad: comprobar que nombres de conteos y metadata coinciden
all(rownames(coldata) == colnames(data))
head(coldata)
str(coldata)

## 3. Visualización inicial previa a la construcción del objeto DESeq2 ####

# Subdata
data_subset <- data %>% 
  dplyr::select(PN1, PN2, PN3, PN4, PN5, PT1, PT2, PT3, PT4, PT5)
dim(data_subset)

data_long <- data_subset %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "counts")  # Transformar a formato largo para que ggplot pueda manejar todas las muestras
head(data_long)

# Histograma dataset no filtrado
hist_plot_pre <- ggplot(data_long, aes(x = counts)) + 
  geom_histogram(bins = 75, fill = "#3C5488", color = "black", linewidth = 0.2) + 
  scale_x_log10() +  
  labs(x = "Raw expression counts (log10)", y = "Number of genes") +
  facet_wrap(~sample, ncol = 5, nrow = 2) +
  tema_paper +
  theme(axis.text.x = element_blank(), strip.background = element_rect(fill="gray90"))
hist_plot_pre

## Guardar Histograma
ggsave(
  filename = file.path(output_dir, "histograma_counts_prefiltered.png"), 
  plot = hist_plot_pre, 
  width = 10, 
  height = 5, 
  dpi = 300
)

# Principal component analysis (PCA)
varianza_genes <- apply(data_subset, 1, var)  # varianza por FILAS (genes)
data_filtered_PCA <- data_subset[varianza_genes > 0.1 & !is.na(varianza_genes),] # Filtrado porque PCA no admite columnas (genes) con todo 0
dim(data_filtered_PCA)

pca.results <- prcomp(t(data_filtered_PCA),       
                      center=TRUE, 
                      scale.=TRUE) 

data.pca <- data.frame(pca.results$x)

data.pca$Sample <- rownames(data.pca)
data.pca$Class <- substr(data.pca$Sample, 1, 2)
data.pca$Patient <- substr(data.pca$Sample, 3, 3)
str(data.pca)

varianzas <- pca.results$sdev^2
varianza.explicada <- varianzas / sum(varianzas)

## PCA by Class
PCAbyclass <- ggplot(data.pca, aes(x=PC1, y=PC2, color=Class)) +  
  geom_point(size=4, alpha=0.8) +  
  scale_color_manual(values=c("PN" = "#00A087", "PT" = "#F39B7F")) + 
  labs(title= "PCA: Normal (PN) vs Tumor (PT)", 
       x=paste0('PC1 (', round(varianza.explicada[1] * 100, 2), '%)'),                    
       y=paste0('PC2 (', round(varianza.explicada[2] * 100, 2),'%)'),                   
       color='Condición') +             
  tema_paper  
PCAbyclass

ggsave(
  filename = file.path(output_dir, "PCA_class_prefiltered.png"), 
  plot = PCAbyclass, 
  width = 10, 
  height = 8, 
  dpi = 300
)

## PCA by Patient
PCAbypatient <- ggplot(data.pca, aes(x=PC1, y=PC2, color=Patient)) +  
  geom_point(size=4, alpha=0.8) +  
  scale_color_brewer(palette = "Set1") + # Paleta apta para papers
  labs(title= "PCA: Variabilidad por Paciente", 
       x=paste0('PC1 (', round(varianza.explicada[1] * 100, 2), '%)'),                    
       y=paste0('PC2 (', round(varianza.explicada[2] * 100, 2),'%)'),                   
       color='Paciente') +             
  tema_paper
PCAbypatient

ggsave(
  filename = file.path(output_dir, "PCA_patient_prefiltered.png"), 
  plot = PCAbypatient, 
  width = 10, 
  height = 8, 
  dpi = 300
)

## 4. Construcción del objeto DESeq2 y Filtrado ####

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ Patient + Condition)
dim(dds) # [1] 58395    15

# Filtrado de genes con baja expresión (al menos 10 counts en el tamaño del grupo más pequeño)

smallestGroupSize <- 5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # [1] 17525    15


## 5. Ejecución de DESeq2 y Normalización ####

# Ejecutar el pipeline principal
dds <- DESeq(dds)

# Transformación logarítmica regularizada (rlog) para Visualización (PCA, Heatmaps)
rld <- rlog(dds, blind=FALSE) # blind=FALSE porque ya conocemos el diseño experimental
rld_mat <- assay(rld)

## 6. Análisis Exploratorio de Datos (EDA) y Control de Calidad ####

# 6.1 Histograma de conteos transformados
data_long_rlog <- as.data.frame(rld_mat) %>%
  dplyr::select(PN1, PN2, PN3, PN4, PN5, PT1, PT2, PT3, PT4, PT5) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "rlog_counts")

hist_plot <- ggplot(data_long_rlog, aes(x = rlog_counts)) + 
  geom_histogram(bins = 75, fill = "#3C5488", color = "black", linewidth = 0.2) + 
  labs(x = "rlog transformed counts", y = "Number of genes") +
  facet_wrap(~sample, ncol = 5, nrow = 3) +
  tema_paper +
  theme(strip.background = element_rect(fill="gray90"))

hist_plot

# Guardar Histograma
ggsave(
  filename = file.path(output_dir, "histograma_counts_filtered.png"), 
  plot = hist_plot, 
  width = 10, 
  height = 5, 
  dpi = 300
)

# 6.2 PCA (Análisis de Componentes Principales)

# Filtrar muestras_PA, solo queremos ver PN vs PT
muestras_sin_PA <- colData(dds)$Condition %in% c("PN", "PT")
rld_sin_PA <- rld[, muestras_sin_PA]
colData(rld_sin_PA)$Condition <- droplevels(colData(rld_sin_PA)$Condition)

# PCA con los top 500 genes más variables
pca_plot500 <- plotPCA(rld_sin_PA, intgroup = c("Condition"), ntop=500) + 
  scale_color_manual(values = c("PN" = "#00A087", "PT" = "#F39B7F"),
                     labels = c("Normal (PN)", "Tumor (PT)")) +
  labs(title = NULL,
       color = "Condición:") +
  tema_paper

pca_plot500


# PCA con los top 1500 genes más variables
pca_plot1500 <- plotPCA(rld_sin_PA, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values = c("PN" = "#00A087", "PT" = "#F39B7F"),
                     labels = c("Normal (PN)", "Tumor (PT)")) +
  labs(title = NULL,
       color = "Condición:") +
  tema_paper

pca_plot1500

# Guardar PCA
ggsave(
  filename = file.path(output_dir, "PCA_filtered.png"), 
  plot = pca_plot1500, 
  width = 10, 
  height = 5, 
  dpi = 300
)

## 7. Análisis de Expresión Diferencial (DGE) ####

?results
resultsNames(dds)

# Extraer resultados crudos (PT vs PN)
res <- results(dds, alpha=0.05, lfcThreshold = 0.585)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE) # 1519 (540 up y 979 down) *mirar final del script

# Encogimiento del Log Fold Change (Shrinkage) para reducir falsos positivos
# y penalizar genes con conteos muy bajos o alta variabilidad
resLFC <- lfcShrink(dds, coef="Condition_PT_vs_PN", type="apeglm")
resLFC

## 7.1. AnnotationHub ####

# Crear el objeto Hub (la conexión)
ah <- AnnotationHub()

# Buscar bases de datos de humano y Ensembl
query(ah, c("Homo sapiens", "EnsDb",
            "103")) # Anotación contemporánea a la publicación original (~ enero de 2021)

# Descargar la base de datos específica
edb <- ah[["AH89426"]]
hubCache(ah)

# Añadir columna ENTREZID a la tabla de resultados
resLFC$entrezid <- mapIds(edb,
                          keys = rownames(resLFC),      # Tus ENSG...
                          column = "ENTREZID",          # Identificador requerido para pathways
                          keytype = "GENEID",           # Lo que tienes en los rownames
                          multiVals = "first")

# Añadir columna SYMBOL a la tabla de resultados
resLFC$symbol <- mapIds(edb,
                        keys = rownames(resLFC),      # Tus ENSG...
                        column = "SYMBOL",            # Lo que quieres obtener
                        keytype = "GENEID",           # Lo que tienes ahora
                        multiVals = "first")          # Si hay dudas, dame el primero
head(resLFC)

#####
# Ordenar por p-valor ajustado
resLFC_ordered <- resLFC[order(resLFC$padj), ]
resLFC_ordered

# Filtrar genes significativos (FoldChange > 1.5x y p-adj < 0.05)
resLFC_filtered <- resLFC_ordered[abs(resLFC_ordered$log2FoldChange) > 0.585 & 
                                    resLFC_ordered$padj < 0.05 & 
                                    !is.na(resLFC_ordered$padj), ]
resLFC_filtered
dim(resLFC_filtered)
sum(is.na(resLFC_filtered)) #9 (symbol) y 459 (entrezid)

# Extraer los Top 20 genes más significativos para visualización posterior
top20_sig_genes <- rownames(resLFC_filtered)[1:20]
top20_sig_symbols <- resLFC_filtered$symbol[1:20]

top1500_sig_genes <- rownames(resLFC_filtered)[1:1500]

# Guardar tabla de resultados diferencialmente expresados
write.csv(as.data.frame(resLFC_filtered), file = file.path(output_dir, "Resultados_LFC_Shrink.csv"))

# Filtrar genes significativos (FoldChange > 2x y p-adj < 0.05)
resLFC_filtered_1 <- resLFC_ordered[abs(resLFC_ordered$log2FoldChange) > 1 & 
                                      resLFC_ordered$padj < 0.05 & 
                                      !is.na(resLFC_ordered$padj), ]

dim(resLFC_filtered_1) # [1] 3789    7


file_name <- paste0("Resultados_LFC_Shrink_", dataset_name, ".csv")
full_path <- file.path("Output/DEGS", file_name)
write.csv(as.data.frame(resLFC_filtered_1), file = full_path)

## 8. Visualización de Resultados de DGE ####

# 8.1 PCA
rld_top1500_sig <- rld_sin_PA[top1500_sig_genes, ]

pca_plot1500_sig <- plotPCA(rld_top1500_sig, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values = c("PN" = "#00A087", "PT" = "#F39B7F"),
                     labels = c("Normal (PN)", "Tumor (PT)")) +
  labs(title = NULL,
       color = "Condición:") +
  tema_paper

pca_plot1500_sig

ggsave(
  filename = file.path(output_dir, "PCA_top1500sign_filtered.png"), 
  plot = pca_plot1500_sig, 
  width = 10, 
  height = 5, 
  dpi = 300
)

# 8.2 MA-Plot

# Guardar MA-Plot usando la interfaz de base R (png/dev.off)
png(file.path(output_dir, "MA_Plot.png"), width = 800, height = 600, res = 120)
par(family = "serif", mar = c(5, 5, 4, 2) + 0.1) # Cambia a Times New Roman y ajusta márgenes
plotMA(resLFC, 
       ylim = c(-6, 6),            
       cex = 0.6,
       colNonSig = "gray60", colSig = "blue", colLine = "black")  
abline(h = c(-0.585, 0.585), col = "#DC0000", lwd = 2, lty = 2) 
dev.off()

# 8.3 Boxplots de los Top 20 Genes
top20_counts <- counts(dds, normalized = TRUE)[top20_sig_genes, ]
top20_df <- as.data.frame(top20_counts) %>%
  mutate(gene = top20_sig_symbols) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "counts") %>%
  left_join(as.data.frame(colData(dds)) %>% mutate(sample = rownames(.)), by = "sample")

# Mantener solo PN y PT
top20_filtered <- top20_df %>% 
  dplyr::filter(Condition %in% c("PN", "PT"))

boxplot_top20 <- ggplot(top20_filtered, aes(x = Condition, y = counts, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, color = "black") +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) + 
  scale_y_log10() + 
  scale_fill_manual(values = c("PN" = "#3C5488", "PT" = "#DC0000")) +
  labs(title = "Top 20 Genes: Normal (PN) vs Tumor (PT)",
       y = "Normalized Counts (log10)",
       x = "Condición Biológica") +
  tema_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "italic", family = "serif"),
        strip.background = element_rect(fill="white", color="black"))
boxplot_top20

ggsave(file.path(output_dir, "Boxplot_Top20_Genes.png"),
       plot = boxplot_top20,
       width = 12, height = 8, dpi = 300)

# 8.4 Heatmap de los Top 20 Genes

rld_sin_PA_mat <- assay(rld_sin_PA)
head(rld_sin_PA_mat)

# Extraer matriz rlog solo de genes top
mat_filtrada <- rld_sin_PA_mat[top20_sig_genes,]
rownames(mat_filtrada) <- top20_sig_symbols
head(mat_filtrada)

# Guardar Heatmap (pheatmap tiene argumento filename integrado)
pheatmap(mat_filtrada, 
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         cluster_cols = TRUE,    
         annotation_col = coldata,
         main = "Top 20 Genes: PN vs PT",
         fontfamily = "serif",     # Asegura la fuente Times New Roman
         fontsize = 10,
         filename = file.path(output_dir, "Heatmap_Top20_Genes.png"),
         width = 8, height = 6)

## 8.5 Volcano Plot ####

# Crear una nueva columna para clasificar la significancia de cada gen
# Umbrales: padj < 0.05 y abs(log2FoldChange) > 0.585 (1.5x)
volcano_data <- as.data.frame(resLFC_ordered) %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 0.585 ~ "Upregulated en PT",
      padj < 0.05 & log2FoldChange < -0.585 ~ "Downregulated en PT",
      TRUE ~ "Not Significant"
    )
  )
head(volcano_data)

# Filtrar los datos para aislar solo los Top 20 genes (para etiquetarlos en el gráfico)
top_genes_volcano <- volcano_data %>% dplyr::filter(symbol %in% top20_sig_symbols)

# Construir el gráfico con ggplot2
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Downregulated en PT" = "#3C5488", 
                                "Not Significant" = "gray80", 
                                "Upregulated en PT" = "#DC0000")) +
  geom_vline(xintercept = c(-0.585, 0.585), col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(data = top_genes_volcano, aes(label = symbol), 
                  size = 3.5, color = "black", family = "serif", # Fuente de las etiquetas
                  fontface = "italic", box.padding = 0.5, max.overlaps = 20) +
  labs(title = "Volcano Plot: Tumor vs Normal",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "") +
  tema_paper

volcano_plot

ggsave(file.path(output_dir, "Volcano_Plot.png"), plot = volcano_plot, width = 8, height = 6, dpi = 300)






## 0. Bloque de Exploración de Umbrales ####

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE) # 4607

res005 <- results(dds, alpha=0.005)
summary(res005)
sum(res005$padj < 0.005, na.rm=TRUE) # 2050

res001 <- results(dds, alpha=0.001)
summary(res001)
sum(res001$padj < 0.001, na.rm=TRUE) # 1238

### Foldchange de x1.5

res05_05 <- results(dds, alpha=0.05, lfcThreshold = 0.585)
summary(res05_05)
sum(res05_05$padj < 0.05, na.rm=TRUE) # 1519

res005_05 <- results(dds, alpha=0.005, lfcThreshold = 0.585)
summary(res005_05)
sum(res005_05$padj < 0.005, na.rm=TRUE) # 665

res001_05 <- results(dds, alpha=0.001, lfcThreshold = 0.585)
summary(res001_05)
sum(res001_05$padj < 0.001, na.rm=TRUE) # 441

### Foldchange de x2

res05_1 <- results(dds, alpha=0.05, lfcThreshold = 1)
summary(res05_1)
sum(res05_1$padj < 0.05, na.rm=TRUE) # 645

res005_1 <- results(dds, alpha=0.005, lfcThreshold = 1)
sum(res005_1$padj < 0.005, na.rm=TRUE) # 327

res001_1 <- results(dds, alpha=0.001, lfcThreshold = 1)
sum(res001_1$padj < 0.001, na.rm=TRUE) # 209


#### Esto elimina todos los archivos descargados por AnnotationHub 
#AnnotationHub::removeCache(ah)
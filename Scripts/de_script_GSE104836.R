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
dataset_name <- "GSE104836"

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
data <- read.table("Data/GSE104836/GSE104836_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)
head(data)
dim(data)

# Carga de metadata
meta <- read.csv("Data/GSE104836/GSE104836_full_metadata.csv", header = TRUE, row.names = 1)
head(meta)

# Generar metadata (coldata) extrayendo información de los nombres de muestra
coldata <- data.frame(
  Condition = as.factor(substr(meta$Sample_characteristics_ch1, 9, nchar(meta$Sample_characteristics_ch1) - 7)),
  Patient = as.factor(substr(meta$Sample_title, 1, 2)),
  row.names = colnames(data)
)

# Relevel para que Normal (nontumor colon tissue) sea la referencia basal
coldata$Condition <- relevel(coldata$Condition, ref = "nontumor colon")

# Check de seguridad: comprobar que nombres de conteos y metadata coinciden
all(rownames(coldata) == colnames(data))
head(coldata)
str(coldata)

## 3. Visualización inicial previa a la construcción del objeto DESeq2 ####

# Transformar a formato largo para que ggplot pueda manejar todas las muestras
data_long <- data %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "counts")  
head(data_long)

# Histograma dataset no filtrado
hist_plot_pre <- ggplot(data_long, aes(x = counts)) + 
  geom_histogram(bins = 75, fill = "#3C5488", color = "black", linewidth = 0.2) + 
  scale_x_log10() +  
  labs(x = "Raw expression counts (log10)", y = "Number of genes") +
  facet_wrap(~sample, ncol = 5, nrow = 4) +
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
varianza_genes <- apply(data, 1, var)  # varianza por FILAS (genes)
data_filtered_PCA <- data[varianza_genes > 0.1 & !is.na(varianza_genes),] # Filtrado porque PCA no admite columnas (genes) con todo 0
dim(data_filtered_PCA)

pca.results <- prcomp(t(data_filtered_PCA),       
                      center=TRUE, 
                      scale.=TRUE) 

data.pca <- data.frame(pca.results$x)

data.pca$Sample <- rownames(data.pca)
data.pca$Class <- coldata$Condition
str(data.pca)

varianzas <- pca.results$sdev^2
varianza.explicada <- varianzas / sum(varianzas)

## PCA by Class
PCAbyclass <- ggplot(data.pca, aes(x=PC1, y=PC2, color=Class)) +  
  geom_point(size=4, alpha=0.8) +  
  scale_color_manual(values=c("nontumor colon" = "#00A087", "colon cancer" = "#F39B7F")) + 
  labs(title= "PCA: Normal (nontumor colon) vs Tumor (colon cancer)", 
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


## 4. Construcción del objeto DESeq2 y Filtrado ####

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ Patient + Condition)
dim(dds) # [1] 39376    20

# Filtrado de genes con baja expresión (al menos 10 counts en el tamaño del grupo más pequeño)
table(coldata$Condition)
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # [1] 23309    20


## 5. Ejecución de DESeq2 y Normalización ####

# Ejecutar el pipeline principal
dds <- DESeq(dds)

# Transformación logarítmica regularizada (rlog) para Visualización (PCA, Heatmaps)
rld <- rlog(dds, blind=FALSE) # blind=FALSE porque ya conocemos el diseño experimental
rld_mat <- assay(rld)

## 6. Análisis Exploratorio de Datos (EDA) y Control de Calidad ####

# 6.1 Histograma de conteos transformados
data_long_rlog <- as.data.frame(rld_mat) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "rlog_counts")

hist_plot <- ggplot(data_long_rlog, aes(x = rlog_counts)) + 
  geom_histogram(bins = 75, fill = "#3C5488", color = "black", linewidth = 0.2) + 
  labs(x = "rlog transformed counts", y = "Number of genes") +
  facet_wrap(~sample, ncol = 5, nrow = 4) +
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

# PCA con los top 500 genes más variables
pca_plot500 <- plotPCA(rld, intgroup = c("Condition"), ntop=500) + 
  scale_color_manual(values=c("nontumor colon" = "#00A087", "colon cancer" = "#F39B7F"),
                     labels = c("Normal (nontumor colon)", "Tumor (colon cancer)")) +
  labs(title = NULL,
       color = "Condición:") +
  tema_paper

pca_plot500


# PCA con los top 1500 genes más variables
pca_plot1500 <- plotPCA(rld, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values=c("nontumor colon" = "#00A087", "colon cancer" = "#F39B7F"),
                     labels = c("Normal (nontumor colon)", "Tumor (colon cancer)")) +
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

# Extraer resultados crudos
res <- results(dds, alpha=0.05, lfcThreshold = 0.585)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE) # 2021 (911 up y 1110 down) *mirar final del script

# Encogimiento del Log Fold Change (Shrinkage) para reducir falsos positivos
# y penalizar genes con conteos muy bajos o alta variabilidad
resLFC <- lfcShrink(dds, coef="Condition_colon.cancer_vs_nontumor.colon", type="apeglm")
resLFC

## 7.1. AnnotationHub ####

# Crear el objeto Hub (la conexión)
ah <- AnnotationHub()

# Buscamos la base de datos OrgDb de humano
busqueda_org <- query(ah, c("Homo sapiens", "OrgDb"))
df_org <- as.data.frame(mcols(busqueda_org))

# Descargamos el objeto usando el ID que encontraste
gdb_humano <- ah[["AH121953"]]

# 3. Creamos la columna 'symbol' en tus resultados (resLFC)
resLFC$symbol <- mapIds(gdb_humano, 
                        keys = rownames(resLFC), 
                        column = "SYMBOL", 
                        keytype = "ENTREZID", # Indica que tus filas son números Entrez
                        multiVals = "first")
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
dim(resLFC_filtered) # [1] 6742    6
sum(is.na(resLFC_filtered)) #216

# Extraer los Top 20 genes más significativos para visualización posterior
top20_sig_genes <- rownames(resLFC_filtered)[1:20]
top20_sig_symbols <- resLFC_filtered$symbol[1:20]

top1500_sig_genes <- rownames(resLFC_filtered)[1:1500]

# Guardar tabla de resultados diferencialmente expresados
write.csv(as.data.frame(resLFC_filtered), file = file.path(output_dir, "Resultados_LFC_Shrink.csv"))

## 8. Visualización de Resultados de DGE ####

# 8.1 PCA
rld_top1500_sig <- rld[top1500_sig_genes, ]

pca_plot1500_sig <- plotPCA(rld_top1500_sig, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values=c("nontumor colon" = "#00A087", "colon cancer" = "#F39B7F"),
                     labels = c("Normal (nontumor colon)", "Tumor (colon cancer)")) +
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

boxplot_top20 <- ggplot(top20_df, aes(x = Condition, y = counts, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8, color = "black") +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) + 
  scale_y_log10() + 
  scale_fill_manual(values=c("nontumor colon" = "#00A087", "colon cancer" = "#F39B7F")) +
  labs(title = "Top 20 Genes: Normal vs Tumor",
       y = "Normalized Counts (log10)",
       x = NULL,            # Eliminamos el título del eje X
       fill = "Condición") + # Título personalizado para la leyenda
  tema_paper +
  theme(
    legend.position = "bottom", # Activamos la leyenda abajo
    strip.text = element_text(size = 10, face = "italic", family = "serif"),
    strip.background = element_rect(fill="white", color="black"),
    # Quitamos las etiquetas y las marcas del eje X
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

boxplot_top20

ggsave(file.path(output_dir, "Boxplot_Top20_Genes.png"),
       plot = boxplot_top20,
       width = 12, height = 8, dpi = 300)

# 8.4 Heatmap de los Top 20 Genes

# Extraer matriz rlog solo de genes top
mat_filtrada <- rld_mat[top20_sig_genes,]
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
sum(res05$padj < 0.05, na.rm=TRUE) # 9325

res005 <- results(dds, alpha=0.005)
summary(res005)
sum(res005$padj < 0.005, na.rm=TRUE) # 5537

res001 <- results(dds, alpha=0.001)
summary(res001)
sum(res001$padj < 0.001, na.rm=TRUE) # 4056

### Foldchange de x1.5

res05_05 <- results(dds, alpha=0.05, lfcThreshold = 0.585)
summary(res05_05)
sum(res05_05$padj < 0.05, na.rm=TRUE) # 2021

res005_05 <- results(dds, alpha=0.005, lfcThreshold = 0.585)
summary(res005_05)
sum(res005_05$padj < 0.005, na.rm=TRUE) # 1166

res001_05 <- results(dds, alpha=0.001, lfcThreshold = 0.585)
summary(res001_05)
sum(res001_05$padj < 0.001, na.rm=TRUE) # 854

### Foldchange de x2

res05_1 <- results(dds, alpha=0.05, lfcThreshold = 1)
summary(res05_1)
sum(res05_1$padj < 0.05, na.rm=TRUE) # 724

res005_1 <- results(dds, alpha=0.005, lfcThreshold = 1)
sum(res005_1$padj < 0.005, na.rm=TRUE) # 429

res001_1 <- results(dds, alpha=0.001, lfcThreshold = 1)
sum(res001_1$padj < 0.001, na.rm=TRUE) # 316


#### Esto elimina todos los archivos descargados por AnnotationHub 
#AnnotationHub::removeCache(ah)
################################################################################
#
#   title: "Gene-level differential expression analysis using DESeq2"
#   author: "Verónica Cabeza de Vaca Tocino"
#   date: "2026-04"
#
################################################################################

## 1. Carga de librerías ####

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(ensembldb)
library(AnnotationDbi)

getwd()

# Crear directorio para guardar resultados si no existe
dataset_name <- "TCGA-COAD"

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

# 2. Importación de los datos y unificación de las matrices de conteo

# Definir las rutas de trabajo
ruta_datos <- "Data/TCGA-COAD" 
ruta_samplesheet <- "Data/TCGA-COAD/gdc_sample_sheet.2026-05-15.tsv"

# Leer la Sample Sheet
meta_raw <- read_tsv(ruta_samplesheet, show_col_types = FALSE)

# Listar todos los archivos .tsv descargados de las subcarpetas
archivos_tsv <- list.files(path = ruta_datos, 
                           pattern = "\\.tsv$", 
                           recursive = TRUE, 
                           full.names = TRUE)

archivos_tsv <- archivos_tsv[!str_detect(archivos_tsv, "sample_sheet")]

# Función para leer un archivo, limpiar estadísticas y extraer los conteos
leer_conteos <- function(ruta_archivo) {
  nombre_archivo <- basename(ruta_archivo)  # Extraer el nombre del archivo
  datos <- read_tsv(ruta_archivo,
                    skip = 1)   # Leer el archivo saltando la primera línea
  # Filtramos y nos quedamos solo con gene_id y unstranded
  conteos <- datos %>%
    dplyr::filter(!str_detect(gene_id, "^N_")) %>%
    dplyr::select(gene_id, unstranded)
  # Renombrar la columna 'unstranded' con el nombre del archivo temporalmente
  colnames(conteos)[2] <- nombre_archivo
  return(conteos)
}

# Aplicar la función a todos los archivos y fusionar (merge) todo por 'gene_id'
lista_conteos <- lapply(archivos_tsv, leer_conteos)
matriz_cruda <- purrr::reduce(lista_conteos, full_join, by = "gene_id")

# Mapear los nombres de los archivos a los IDs de los pacientes (Sample IDs)
nombres_archivos <- colnames(matriz_cruda)[-1] # Descartar la primera columna (gene_id)

# Crear un vector diccionario: File Name -> Sample ID
diccionario_nombres <- meta_raw$`Sample ID`
names(diccionario_nombres) <- meta_raw$`File Name`

# Aplicar el cambio de nombres a la matriz
nuevos_nombres <- diccionario_nombres[nombres_archivos]
colnames(matriz_cruda)[-1] <- nuevos_nombres

# Formatear como matriz estándar de R 
matriz_final <- as.data.frame(matriz_cruda)
rownames(matriz_final) <- matriz_final$gene_id # Poner los genes como nombre de fila
matriz_final$gene_id <- NULL                   # Borrar la columna gene_id
data <- as.data.frame(matriz_final)        # Convertir a formato matriz

# Echar un vistazo a las primeras 5 filas y 5 columnas para comprobar
head(data[, 1:5])

# Eliminar duplicados del 'Sample ID' 
meta <- meta_raw %>%
  distinct(`Sample ID`, .keep_all = TRUE)

# Convertir a data.frame  y asignar los rownames de forma segura
meta <- as.data.frame(meta)
rownames(meta) <- meta$`Sample ID`
head(meta)

# Filtramos y ordenamos las columnas de la matriz para que coincidan EXACTAMENTE con los rownames de meta
data <- data[, rownames(meta)]
dim(data) # 60660   514

# Comprobación de seguridad para DESeq2 (debería devolver TRUE)
all(colnames(data) == rownames(meta))

## 3. Preparación de Datos para DESeq2 ####

# Generar coldata extrayendo información de los nombres de muestra
coldata <- data.frame(
  Condition = as.factor(meta$'Tissue Type'),
  row.names = rownames(meta)
)

# Relevel para que Normal (Normal) sea la referencia basal
coldata$Condition <- relevel(coldata$Condition, ref = "Normal")
table(coldata$Condition)

# Check de seguridad: comprobar que nombres de conteos y metadata coinciden
all(rownames(coldata) == colnames(data))
head(coldata)
str(coldata)

## 4. Visualización inicial previa a la construcción del objeto DESeq2 ####


# Subdata
data_long <- data %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "counts")  

coldata_limpio <- coldata %>%
  rownames_to_column(var = "sample")

data_long_final <- data_long %>%
  left_join(coldata_limpio, by = "sample")

head(data_long_final)

# Gráfico de densidad agrupando por Condición
density_plot <- ggplot(data_long_final, aes(x = log2(counts + 1), group = sample, color = Condition)) +
  geom_density(alpha = 0.2, linewidth = 0.2) + 
  scale_color_manual(values = c("Normal" = "#00A087", "Tumor" = "#F39B7F")) +
  labs(x = "Raw expression counts (Log2)", y = "Density", title = "Distribución de lecturas por muestra") +
  tema_paper 

density_plot

ggsave(
  filename = file.path(output_dir, "densityplot_counts_prefiltered.png"), 
  plot = density_plot, 
  width = 8, 
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
  geom_point(size=4, alpha=0.5) +  
  scale_color_manual(values=c("Normal" = "#00A087", "Tumor" = "#F39B7F")) + 
  labs(title= "PCA: Normal vs Tumor", 
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
                              design = ~ Condition)
dim(dds) # [1] 39376    54

# Filtrado de genes con baja expresión (al menos 10 counts en el tamaño del grupo más pequeño)
table(coldata$Condition)
smallestGroupSize <- 41
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dim(dds) # [1] 24160   514

genes_universo_enrichGO <- rownames(dds)
  
## 5. Ejecución de DESeq2 y Normalización ####

# Ejecutar el pipeline principal
dds <- DESeq(dds)

# Transformación logarítmica regularizada (rlog) para Visualización (PCA, Heatmaps)
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)

## 6. Análisis Exploratorio de Datos (EDA) y Control de Calidad ####

# 6.1 Gráfico de densidad agrupando por Condición
data_long_vst <- as.data.frame(vsd_mat) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "vst_counts") %>%
  left_join(as.data.frame(colData(dds)) %>% rownames_to_column(var = "sample"), by = "sample")

density_plot_vst <- ggplot(data_long_vst, aes(x = vst_counts, group = sample, color = Condition)) +
  geom_density(alpha = 0.1, linewidth = 0.2) + # Bajamos un poco el grosor por ser 160 líneas
  scale_color_manual(values = c("Normal" = "#00A087", "Tumor" = "#F39B7F")) +
  labs(x = "VST Transformed Counts", 
       y = "Density", 
       title = "Distribución de lecturas normalizadas") +
  tema_paper 

density_plot_vst

ggsave(
  filename = file.path(output_dir, "densityplot_counts_filtered.png"), 
  plot = density_plot_vst, 
  width = 8, 
  height = 5, 
  dpi = 300
)

# 6.3 PCA (Análisis de Componentes Principales)

# PCA con los top 500 genes más variables
pca_plot500 <- plotPCA(vsd, intgroup = c("Condition"), ntop=500) + 
  scale_color_manual(values = c("Normal" = "#00A087", "Tumor" = "#F39B7F"),
                     labels = c("Normal", "Tumor")) +
  labs(title = NULL,
       color = "Condición:") +
  tema_paper

pca_plot500


# PCA con los top 1500 genes más variables
pca_plot1500 <- plotPCA(vsd, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values = c("Normal" = "#00A087", "Tumor" = "#F39B7F"),
                     labels = c("Normal", "Tumor")) +
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
sum(res$padj < 0.05, na.rm=TRUE) # 8364 (5021 up y 3343 down) *mirar final del script

# Encogimiento del Log Fold Change (Shrinkage) para reducir falsos positivos
# y penalizar genes con conteos muy bajos o alta variabilidad
resLFC <- lfcShrink(dds, coef="Condition_Tumor_vs_Normal", type="apeglm", res=res)
resLFC

## 7.1. Inclusión de nomenclatura symbol: Adoptar gene symbol de .tsv originales ####

archivo_diccionario <- archivos_tsv[1]

# Leemos el archivo saltando la primera línea de metadatos, 
# filtramos las filas de estadísticas y nos quedamos con ambas columnas identificadoras.
diccionario_genes <- read_tsv(archivo_diccionario, show_col_types = FALSE, skip = 1) %>%
  dplyr::filter(!str_detect(gene_id, "^N_")) %>%
  dplyr::select(gene_id, symbol = gene_name) # Renombramos aquí + Aseguramos que no haya IDs de Ensembl duplicados 

# Preparar el objeto de resultados de DESeq2 (resLFC)
resLFC_df <- as.data.frame(resLFC) %>%
  # Pasamos los rownames a una columna real para poder hacer el cruce
  rownames_to_column(var = "gene_id") 

# Cruzar los datos (Matching)
resLFC_df <- resLFC_df %>%
  right_join(diccionario_genes, by = "gene_id")

# Comprobar cómo ha quedado la tabla final
head(resLFC_df)


## 7.2. AnnotationHub para TCGA-COAD ####

ah <- AnnotationHub()

# Buscamos la base de datos de humano. 
busqueda_org_tcga <- query(ah, c("Homo sapiens", "EnsDb"))
df_org_tcga <- as.data.frame(mcols(busqueda_org_tcga))
View(df_org_tcga)

# Descargamos el objeto de Ensembl 103 (equivalente a GENCODE v36 para TCGA-COAD)
gdb_tcga <- ah[["AH89426"]]

# Inclusión de nomenclatura entrezid ####
# PASO CLAVE: Eliminar la versión (GENCODE v36) del ID de Ensembl -> Esto transforma "ENSG00000141510.17" en "ENSG00000141510" para poder cruzarlo
resLFC_df <- resLFC_df %>%
  mutate(gene_id_limpio = str_remove(gene_id, "\\..*$"))

# Mapear de Ensembl a Entrez ID utilizando la base de datos de AnnotationHub específica para TCGA
resLFC_df$entrezid <- mapIds(
  x = gdb_tcga,              # Usamos el objeto descargado para TCGA-COAD en lugar del estático org.Hs.eg.db
  keys = resLFC_df$gene_id_limpio,
  column = "ENTREZID",       # Lo que queremos obtener
  keytype = "GENEID",       # Lo que le pasamos
  multiVals = "first"        # Si hay colisiones (1 a varios), nos quedamos con el primero
)

rownames(resLFC_df) <- resLFC_df$gene_id
resLFC_df <- dplyr::select(resLFC_df, c(baseMean, log2FoldChange, lfcSE, pvalue, padj, entrezid, symbol))

head(resLFC_df)

# Creación universo para Enrich GO

genes_universo_enrichGO <- as.data.frame(genes_universo_enrichGO) %>%
  dplyr::rename(gene_id = 1) %>%
  mutate(gene_id = str_remove(gene_id, "\\..*$"))

head(genes_universo_enrichGO)

genes_universo_enrichGO$entrezid <- mapIds(
  x = org.Hs.eg.db,
  keys = genes_universo_enrichGO$gene_id,
  column = "ENTREZID",     # Lo que queremos obtener
  keytype = "ENSEMBL",     # Lo que le pasamos
  multiVals = "first"      # Si hay colisiones (1 a varios), nos quedamos con el primero
)

head(genes_universo_enrichGO)
sum(!is.na(genes_universo_enrichGO$entrezid)) #19541
#Al final no se usa esta lista como universo porque la traducción implica una disminucion de las dimensiones
#Por lo que se opta por usar el dataset GSE104836, que aporta un total de 23309 genes + no necesitar traducción 
#write.csv(genes_universo_enrichGO, "Output/Interpretación/universo_genes_enrichGO.csv", row.names = FALSE)

#####
# Ordenar por p-valor ajustado
resLFC_ordered <- resLFC_df[order(resLFC_df$padj), ]
resLFC_ordered

# Filtrar genes significativos (FoldChange > 1.5x y p-adj < 0.05)
resLFC_filtered <- resLFC_ordered[abs(resLFC_ordered$log2FoldChange) > 0.585 & 
                                    resLFC_ordered$padj < 0.05 & 
                                    !is.na(resLFC_ordered$padj), ]
resLFC_filtered
dim(resLFC_filtered) # [1] 8341    7
sum(is.na(resLFC_filtered$entrezid)) #2398
sum(is.na(resLFC_filtered$symbol)) #0

# Extraer los Top 20 genes más significativos para visualización posterior
top20_sig_genes <- rownames(resLFC_filtered)[1:20]
top20_sig_symbols <- resLFC_filtered$symbol[1:20]

top1500_sig_genes <- rownames(resLFC_filtered)[1:1500]

# Guardar tabla de resultados diferencialmente expresados
write.csv(as.data.frame(resLFC_filtered), file = file.path(output_dir, "Resultados_LFC_Shrink.csv"))

## 8. Visualización de Resultados de DGE ####

# 8.1 PCA
vsd_top1500_sig <- vsd[top1500_sig_genes, ]

pca_plot1500_sig <- plotPCA(vsd_top1500_sig, intgroup = c("Condition"), ntop=1500) + 
  scale_color_manual(values=c("Normal" = "#00A087", 
                              "Tumor" = "#F39B7F")) +
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
  facet_wrap(~ gene, scales = "free_y", ncol = 5) + 
  scale_y_log10() + 
  scale_fill_manual(values=c("Normal" = "#00A087", "Tumor" = "#F39B7F")) +
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

mat_filtrada <- assay(vsd[top20_sig_genes, ]) 
rownames(mat_filtrada) <- top20_sig_symbols
head(mat_filtrada)

anno_select <- coldata[, "Condition", drop = FALSE]

pheatmap(mat_filtrada, 
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         cluster_cols = TRUE,    
         annotation_col = anno_select, # Usamos solo la columna Condition
         main = "Top 20 Genes: Normal vs Tumor",
         fontfamily = "serif",
         fontsize = 13,
         fontsize_col = 6,           
         show_colnames = FALSE, #
         filename = file.path(output_dir, "Heatmap_Top20_Genes.png"),
         width = 14, height = 8) 

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
sum(res05$padj < 0.05, na.rm=TRUE) # 17949

res005 <- results(dds, alpha=0.005)
summary(res005)
sum(res005$padj < 0.005, na.rm=TRUE) # 15464

res001 <- results(dds, alpha=0.001)
summary(res001)
sum(res001$padj < 0.001, na.rm=TRUE) # 14090

### Foldchange de x1.5

res05_05 <- results(dds, alpha=0.05, lfcThreshold = 0.585)
summary(res05_05)
sum(res05_05$padj < 0.05, na.rm=TRUE) # 8364

res005_05 <- results(dds, alpha=0.005, lfcThreshold = 0.585)
summary(res005_05)
sum(res005_05$padj < 0.005, na.rm=TRUE) # 6882

res001_05 <- results(dds, alpha=0.001, lfcThreshold = 0.585)
summary(res001_05)
sum(res001_05$padj < 0.001, na.rm=TRUE) # 6123

### Foldchange de x2

res05_1 <- results(dds, alpha=0.05, lfcThreshold = 1)
summary(res05_1)
sum(res05_1$padj < 0.05, na.rm=TRUE) # 4744

res005_1 <- results(dds, alpha=0.005, lfcThreshold = 1)
sum(res005_1$padj < 0.005, na.rm=TRUE) # 3862

res001_1 <- results(dds, alpha=0.001, lfcThreshold = 1)
sum(res001_1$padj < 0.001, na.rm=TRUE) # 3379


#### Esto elimina todos los archivos descargados por AnnotationHub 
#AnnotationHub::removeCache(ah)
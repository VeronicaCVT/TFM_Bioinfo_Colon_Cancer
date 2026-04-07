
setwd("C:/Users/verit/Documents/TFM_Data") # NO debe ser muy alrgo sino no se descarga: file not found


# Step 1: Install and Load TCGAbiolinks

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(Biobase)

dir.create("GDCdata", recursive = TRUE)

# Step 2: View Available TCGA Projects

projects <- getGDCprojects()
head(projects$project_id)


# Step 3: Query the Data

query <- GDCquery(
  project = "TCGA-COAD", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  access = "open" 
)


# Step 4: Download the Data

GDCdownload(query, method = "api", files.per.chunk = 10)


# Step 5: Prepare the Data

data_se <- GDCprepare(query)


# Step 6: Explore and Extract the Count Matrix

head(colData(data_se))

count_matrix <- assay(data_se, "unstranded")  # "unstranded" es la estandar para DESeq2
head(count_matrix)

## To get clinical data:

metadata_final <- as.data.frame(colData(data_se))
metadata_plano <- as.data.frame(metadata_final)

for (col in colnames(metadata_plano)) {
  if (is.list(metadata_plano[[col]])) {
    metadata_plano[[col]] <- sapply(metadata_plano[[col]], paste, collapse = ";")
  }
}

nrow(metadata_final)

clinical_data <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
head(clinical_data)


# DOWNLOAD

# 1. Guardar la matriz de conteos crudos (Genes en filas, Muestras en columnas)
# Usamos row.names = TRUE para no perder los IDs de los genes (Ensembl)
write.csv(count_matrix, file = "TCGA-COAD_Raw_Counts.csv", row.names = TRUE, na = "NA")

# 2. Guardar los metadatos clínicos (Información de los pacientes)
# Convertimos el objeto colData a un data.frame normal para que sea legible
write.csv(metadata_plano, file = "TCGA-COAD_METADATA_Limpio.csv", row.names = FALSE)
write.csv(clinical_data, file = "TCGA-COAD_Clinical_Data.csv", row.names = FALSE, na = "NA")


# DESCRIPCIÓN: 
# Seleccionar las muestras que tienen metilación y transcriptoma

# Settings ----------------------------------------------------------------

rm(list=ls())

# Working directory
root = "/data/local/mvaquerizo/transcriptomica/"
setwd(root)

# Librerías
library(readxl)
library(ggvenn)

# Funciones 
source("Funciones_R.R")

# Carpeta de resultados
results_folder = paste0(root, "results/05_select_transcriptome_samples/")
dir.create(results_folder)

# Carga de datos ---------------------------------------------------------------
# Cargamos el excel de muestras secuenciadas
muestras_secuenciadas = data.frame(read_xlsx(paste0(root, "data/proyectos/Muestras_secuenciadas.xlsx"),
                                             col_names = TRUE))
str(muestras_secuenciadas)
rownames(muestras_secuenciadas) = paste(muestras_secuenciadas$PROJECT, muestras_secuenciadas$ID, sep = "_")

# Cargamos los datos de metilación 
load("results/02_exploratory_analysis/data_filtered.RData", verbose = T)
# Loading objects:
#   Betas
#   ed_exp
#   ed_pat
load("results/04_reproductive_outcome/datos_E_R_Normal.RData", verbose = T)
# Loading objects:
#   sel_datos

# Cargamos datos del TED 
load("../metilacion/results/04_Beta_TED/TED_RUNS_1_2_3_4.RData", verbose = T)
# Loading objects:
#   df_final

# Intersección metilación y transcripción --------------------------------------

# Muestras con transcriptoma
nrow(muestras_secuenciadas[muestras_secuenciadas$TRANSCRIPTOME == 1, ])
# 340 
# Muestras con metilación y transcriptoma
nrow(muestras_secuenciadas[muestras_secuenciadas$TRANSCRIPTOME == 1 & muestras_secuenciadas$METHYLATION == 1, ])
x <- list("Transcriptoma" = rownames(muestras_secuenciadas)[muestras_secuenciadas$TRANSCRIPTOME == 1], 
          "Metilación" = rownames(muestras_secuenciadas)[muestras_secuenciadas$METHYLATION == 1])
ggvenn(x, fill_color = c("#66CDAA", "#CD9ABC"), stroke_size = 0.5,
       set_name_size = 8, text_size = 10) 
# 151 mujeres con ambos datos 

# Seleccionamos las muestras con amba información 
sel = muestras_secuenciadas[muestras_secuenciadas$TRANSCRIPTOME == 1 & muestras_secuenciadas$METHYLATION == 1, ]
muestras_completas = ed_pat[(rownames(ed_pat)%in%sel$Nombre.en.metilación), ]
dim(muestras_completas)
# 148  35
sum(muestras_completas$STUDY == "REPROCYCLE")
# 19 
sum(muestras_completas$STUDY %in% c("TED", "PIF"))
# 77
sum(muestras_completas$STUDY == "ENCLIVA")
# 52
sum(muestras_completas$STUDY == "ENDOAGE")
# 1


# Aplicando los mismos filtros de exclusión de factor embrionario que en metilación
x <- list("Transcriptoma + metilación" = rownames(muestras_completas)[muestras_completas$STUDY != "REPROCYCLE"], "Metilación filtrado" = rownames(ed2_pat_pat))
ggvenn(x, fill_color = c("#66CDAA", "#CD9ABC"), stroke_size = 0.5,
       set_name_size = 8, text_size = 10) 
# 19 controles y 101 pacientes


muestras_transcriptoma = muestras_completas$identificador
muestras_transcriptoma_filtrado = intersect(rownames(muestras_completas), rownames(sel_datos))

filt_muestras_pat = muestras_completas[muestras_transcriptoma_filtrado, ]

# Combinar archivos conteos -----------------------------------------------

# Leemos los ficheros con los counts 
counts_RC1 = read.delim(paste0(root, "data/counts/counts_final_RC_RUN1.tsv"), 
                        sep = "\t", header = T)
rownames(counts_RC1) = gsub("_.*", "", counts_RC1$Geneid)
counts_RC2 = read.delim(paste0(root, "data/counts/counts_final_RC_RUN2.tsv"), 
                        sep = "\t", header = T)
rownames(counts_RC2) = gsub("_.*", "", counts_RC2$Geneid)
counts_PIF = read.delim(paste0(root, "data/counts/gene_counts_noctrl_filtered_PIF.csv"), 
                        sep = ",", header = T)
rownames(counts_PIF) = counts_PIF$X
counts_ENC1 = read.delim(paste0(root, "data/counts/counts_final_ENC_RUN1.tsv"), 
                         sep = "\t", header = T)
rownames(counts_ENC1) = gsub("_.*", "", counts_ENC1$Geneid)
counts_ENC2 = read.delim(paste0(root, "data/counts/counts_final_ENC_RUN2.tsv"), 
                         sep = "\t", header = T)
rownames(counts_ENC2) = gsub("_.*", "", counts_ENC2$Geneid)
counts_ENC3 = read.delim(paste0(root, "data/counts/counts_final_ENC_RUN3.tsv"), 
                         sep = "\t", header = T)
rownames(counts_ENC3) = gsub("_.*", "", counts_ENC3$Geneid)

# Combinamos los arhivos del mismo proyecto 
counts_RC = cbind(counts_RC1, counts_RC2)
counts_ENC = cbind(counts_ENC1, counts_ENC2, counts_ENC3)

#### PIF ####
# Seleccionamos las muestras de PIF --> 57
muestras_PIF = filt_muestras_pat[filt_muestras_pat$STUDY%in%c("TED", "PIF"), ]
ide = muestras_PIF$ID

# Cambiamos nombres columnas transcriptoma, eliminando el run 
colnames(counts_PIF) = gsub("_.*", "", colnames(counts_PIF))

x <- list("PIF" = colnames(counts_PIF), "Metilación" = ide)
ggvenn(x, fill_color = c("skyblue", "seagreen"), stroke_size = 0.5,
       set_name_size = 8, text_size = 10) 

ide[!(ide%in%colnames(counts_PIF))]
# "V55"  "M31"  "P103" "P107" "P98"  "P59"  "V69"  "B12"
# Hay algunas muestras que no están aquí secuenciadas, pero están en los datos de encliva 

# Nos quedamos solo con las columnas comunes con los datos de metilación 
sel_counts_PIF = counts_PIF[, colnames(counts_PIF)%in%ide]
colnames(sel_counts_PIF) = paste("PIF", colnames(sel_counts_PIF), sep = "_")
dim(sel_counts_PIF)
# 20802    47


#### ENCLIVA ####
# Seleccionamos muestras Encliva --> 43
muestras_ENC = filt_muestras_pat[filt_muestras_pat$STUDY == "ENCLIVA", ]
ide = muestras_ENC$identificador

colnames(counts_ENC) = gsub("\\.", "_", colnames(counts_ENC))

x <- list("Encliva" = colnames(counts_ENC), "Metilación" = ide)
ggvenn(x, fill_color = c("skyblue", "seagreen"), stroke_size = 0.5,
       set_name_size = 8, text_size = 10) 
ide[!(ide%in%colnames(counts_ENC))]
# ENC_R9 --> "EMC_R9"
# ENC_Bi34 --> "ENC_Bi34_P1"
# ENC_R49 --> "ENC_R49_P2"

# PIF 
# PIF_V55 --> "PIF_V55"
# PIF_M31 --> "PIF_M31_V2_RNA"
# PIF_P103 --> "PIF_P103_RNA"
# PIF_P107 --> "PIF_P107_RNA"
# PIF_P98 --> "PIF_P98"        "PIF_P98_V2_RNA"
# PIF_P59 --> "PIFP59_V2_RNA"
# PIF_V69 --> "PIF_V69"
# PIF_B12 --> "PIF_B12_V2_RNA"

# EA_V1 --> "EA_V1"

sel_counts_ENC = counts_ENC[, colnames(counts_ENC)%in%ide]

# Añadimos las que faltan 
add = c("EMC_R9", "ENC_Bi34_P1", "ENC_R49_P2", 
        "PIF_V55", "PIF_M31_V2_RNA", "PIF_P103_RNA", "PIF_P107_RNA", 
        "PIF_P98_V2_RNA", "PIFP59_V2_RNA", "PIF_V69", "PIF_B12_V2_RNA", 
        "EA_V1")
sel2_counts_ENC = counts_ENC[, colnames(counts_ENC)%in%add]
colnames(sel2_counts_ENC) = c("ENC_R9", "ENC_Bi34", "ENC_R49",
                              "PIF_V55", "PIF_M31", "PIF_P103", "PIF_P107", 
                              "PIF_P98", "PIF_P59", "PIF_V69", "PIF_B12",
                              "EA_V1")


#### Controles ####
# De controles hay varias porciones de cada paciente
# Para seleccionar cargamos info experimental 
experimental_RC = read.csv(paste0(root, "data/proyectos/Experimental_ReproCycle.csv"), header = T)

# Quitamos controles 
rc_c1 = experimental_RC[grep("RC", experimental_RC$ID), ]
# Nos quedamos solo con el ciclo 1
rc_c1 = rc_c1[grep("C1", rc_c1$ID), ]
# La RC 18 es de ciclo 2 asi que eliminamos el 18 C1 y luego añadimos 18 C2
rc_c1 = rc_c1[grep("18", rc_c1$ID, invert = T), ]
add = experimental_RC[grep("18-C2", experimental_RC$ID), ]
# Eliminamos la RC 13
rc_c1 = rc_c1[grep("13", rc_c1$ID, invert = T), ]
rc = rbind(rc_c1, add)

# Nos quedamos con las muestras que mayor calidad RIN tengan 
sel_rc = rc %>%
  group_by(PACIENTE) %>%
  slice_max(order_by = DV200, n = 1, with_ties = F)

# Unificamos identificadores
sel_rc$ID = gsub("-", "_", sel_rc$ID)
colnames(counts_RC) = gsub("\\.", "_", colnames(counts_RC)) 
all(sel_rc$ID%in%colnames(counts_RC))
# TRUE

# Seleccionamos las columnas del count
sel_counts_RC = counts_RC[, sel_rc$ID]
muestras_RC = colnames(sel_counts_RC)
colnames(sel_counts_RC) = gsub("_C1.*", "", colnames(sel_counts_RC))
colnames(sel_counts_RC) = gsub("_C2.*", "", colnames(sel_counts_RC))

#### Unir ####  
final_counts = data.frame(GeneID = counts_ENC$Geneid)
final_counts =  cbind(final_counts, sel_counts_PIF, sel2_counts_ENC, sel_counts_ENC, sel_counts_RC)
dim(final_counts)
# 20802 119
# 99 pacientes, 19 controles y una columna con el identificador largo 

pat_final_counts = cbind(sel_counts_PIF, sel2_counts_ENC, sel_counts_ENC)
control_final_counts = sel_counts_RC


# Metadata ----------------------------------------------------------------

sel_ed_pat = ed_pat[ed_pat$identificador %in% colnames(pat_final_counts), ]

metadata = data.frame("id_met" = rownames(sel_ed_pat))
metadata$id_trans = sel_ed_pat$identificador
metadata$age = sel_ed_pat$AGE
metadata$bmi = sel_ed_pat$BMI
sel_ed_pat$TR_OUTCOME[sel_ed_pat$TR_OUTCOME %in% c("CM", "BM")] = "M"
metadata$outcome = factor(sel_ed_pat$TR_OUTCOME, levels = c("LB", "M", "NB"))
metadata$p4_hours = sel_ed_pat$P4_HOURS

# Experimental 
sel_ed_exp = ed_exp[metadata$id_met, ]
all(rownames(sel_ed_exp) == metadata$id_met)
# TRUE 
metadata$DV200 = sel_ed_exp$RNA_DV200
metadata$RIN = sel_ed_exp$RIN

# Proyecto 
metadata$project = factor(ifelse(startsWith(metadata$id_trans, "PIF"), "PIF", "ENC"), levels = c("PIF", "ENC"))

# Paciente o control 
metadata$g_paciente = factor(ifelse(startsWith(metadata$id_trans, "RC"), "control", "paciente"),
                             levels = c("paciente", "control"))

# Envejecidas vs rejuvenecidas 
rownames(metadata) = metadata$id_met
sub_sel_datos = sel_datos[metadata$id_met, ]
metadata$g_acc = factor(sub_sel_datos$g_acc, 
                         levels = c("E", "R", "Normal"))


# Las horas de progesterona de PIF_V40 están mal, lo corrijo a mano 
rownames(metadata) = metadata$id_trans
metadata["PIF_V40", "p4_hours"] = difftime("2019-10-10 11:30:00", "2019-10-05 20:00:00", units = "hours")
# Cambiamos el formato a número 
metadata$p4_hours = as.numeric(metadata$p4_hours)

# Otras variables
metadata$estado_embrion  = sel_ed_pat$TR_EMBRYO_STATE
metadata$calidad_embrion = sel_ed_pat$TR_EMBRYO_QUALITY
metadata$origen_ovocito = sel_ed_pat$TR_OOCYTE_ORIGIN
metadata$estado_ovocito = sel_ed_pat$TR_OOCYTE_STATE
metadata$pgt = sel_ed_pat$TR_PGT
metadata$pub = sel_ed_pat$BUP
metadata$ciclo = sel_ed_pat$TR_CYCLE_TYPE

metadata_pat = metadata[metadata$g_acc != "Normal", ]
pat_final_counts = pat_final_counts[, metadata_pat$id_trans]




# Save --------------------------------------------------------------------
save(pat_final_counts, metadata_pat, file = paste0(results_folder, "conteos_finales_E_R.RData"))



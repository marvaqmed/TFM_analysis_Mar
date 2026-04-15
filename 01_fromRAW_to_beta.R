# DESCRIPCIÓN: OBTENCIÓN DE DATOS BETA A PARTIR DE LOS IDATS.
# Cargamos los datos del chip de metilación obtenidos en IIS La Fe. Se generó 
# un data set: "ed".
# Creamos el objeto "methylset" a partir de los IDATs para obtener los valores 
# de los canales RED y GREEN
# Se midieron los parámetros de calidad con el paquete "minfi" e "Ewastools".
# Se obtuvieron los valores Beta a partir de los IDATs, se realizó la media de 
# los sitios Cpgs duplicados para quedarnos solo con la media de ellos. 
# Se cargaron los CSVs de los datos clínicos y experimentales
# Se halló la edad continua AGE2

# Settings ---------------------------------------------------------------------

rm(list=ls()) 

# Working directory
data_folder = "/data/network/nas03/bioinfo_data/fivibio-data/projects/ENDOAGE/"
root = "/data/local/mvaquerizo/metilacion/"
setwd(root)


# Librerías
library(minfi)
library(ewastools)
library(gridExtra)
library(ggplot2)
library(BiocManager)
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
library("IlluminaHumanMethylationEPICv2manifest")
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")

# Funciones 
source("../TFM/Funciones.R")

# Carpeta de resultados
results_folder = paste0(root, "results/01_fromRAW_to_beta/")
dir.create(results_folder)


# Leer sample sheet ------------------------------------------------------------

# Cargamos los datos del chip de metilación de las muestras de cada RUN
# RUN1
samps1 = read.delim(paste0(data_folder, "data/RUN1/RUN1_sampleSheet.csv"), 
                    sep = ",", as.is = T, skip = 7)
# RUN2
samps2 = read.delim(paste0(data_folder, "data/RUN2/RUN2_sampleSheet.csv"), 
                    sep = ",", as.is = T, skip = 7)
#RUN 3
samps3 = read.delim(paste0(data_folder, "data/RUN3/RUN3_sampleSheet.csv"), 
                    sep = ",", as.is = T, skip = 7)
# #RUN 4
samps4 = read.delim(paste0(data_folder, "data/RUN4/RUN4_sampleSheet.csv"), 
                    sep = ",", as.is = T, skip = 7)


# Se genera el data frame con los datos experimentales de la metilación,  
# añadimos el ID de las muestras.
# RUN1
ed1 = data.frame(Sample_Name = paste("PAT", samps1$Sample_Name, sep = "_"),
                Sample_Well = samps1$Sample_Well,
                Sample_Plate = samps1$Sample_Plate,
                Sample_Group = samps1$Type,
                Pool_ID = samps1$Pool_ID,
                Array = samps1$Sentrix_Position,
                Slide = samps1$Sentrix_ID,
                Basename = paste0(data_folder, "data/RUN1/IDATs/", paste(samps1$Sentrix_ID, samps1$Sentrix_Position, sep = "_")), 
                id_corres = paste(samps1$Sentrix_ID, samps1$Sentrix_Position, sep = "_"),
                stringsAsFactors = F)

# RUN2 
ed2 = data.frame(Sample_Name = paste("PAT", samps2$Sample_Name, sep = "_"),
                 Sample_Well = samps2$Sample_Well,
                 Sample_Plate = samps2$Sample_Plate,
                 Sample_Group = samps2$Type,
                 Pool_ID = samps2$Pool_ID,
                 Array = samps2$Sentrix_Position,
                 Slide = samps2$Sentrix_ID,
                 Basename = paste0(data_folder, "data/RUN2/IDATs/", paste(samps2$Sentrix_ID, samps2$Sentrix_Position, sep = "_")), 
                 id_corres = paste(samps2$Sentrix_ID, samps2$Sentrix_Position, sep = "_"),
                 stringsAsFactors = F)

# RUN3 
ed3 = data.frame(Sample_Name = paste("PAT", samps3$Sample_Name, sep = "_"),
                 Sample_Well = samps3$Sample_Well,
                 Sample_Plate = samps3$Sample_Plate,
                 Sample_Group = samps3$Type,
                 Pool_ID = samps3$Pool_ID,
                 Array = samps3$Sentrix_Position,
                 Slide = samps3$Sentrix_ID,
                 Basename = paste0(data_folder, "data/RUN3/IDATs/", paste(samps3$Sentrix_ID, samps3$Sentrix_Position, sep = "_")), 
                 id_corres = paste(samps3$Sentrix_ID, samps3$Sentrix_Position, sep = "_"),
                 stringsAsFactors = F)

# RUN4 
ed4 = data.frame(Sample_Name = paste("PAT", samps4$Sample_Name, sep = "_"),
                 Sample_Well = samps4$Sample_Well,
                 Sample_Plate = samps4$Sample_Plate,
                 Sample_Group = samps4$Type,
                 Pool_ID = samps4$Pool_ID,
                 Array = samps4$Sentrix_Position,
                 Slide = samps4$Sentrix_ID,
                 Basename = paste0(data_folder, "data/RUN4/IDATs/", paste(samps4$Sentrix_ID, samps4$Sentrix_Position, sep = "_")), 
                 id_corres = paste(samps4$Sentrix_ID, samps4$Sentrix_Position, sep = "_"),
                 stringsAsFactors = F)

# Unificar ed de los 4 runs 
ed = rbind(ed1, ed2, ed3, ed4)
head(ed)

# Aquí nombramos los controles: metilado "m", no metilado "u".
ed[ed$Sample_Group=="CTRL_u", "Sample_Name"] = "CTRL_u"
ed[ed$Sample_Group=="CTRL_m", "Sample_Name"] = "CTRL_m"

rownames(ed) = paste(ed$Slide, ed$Array, sep = "_")


# Crear objeto methylSet -------------------------------------------------------

dat = read.metharray.exp(targets = ed)
dat@annotation <- c(array = "IlluminaHumanMethylationEPICv2", 
                    annotation = "20a1.hg38")
dat2 = preprocessRaw(dat)


# Control de calidad -----------------------------------------------------------

# QC minfi : ratio methylated/no_methylated de cada muestra    
dat_QC = minfiQC(dat2, fixOutliers = FALSE, verbose = TRUE)
qc <- as.data.frame(getQC(dat2))
plotQC(qc)

# QC minfi: Densidad de los valores Beta de cada muestra
densityPlot(dat2, sampGroups = ed$Sample_Group)

# En el RUN4  hay una muestra que se va, vamos a intentar ver cual es
ed$fila = gsub("C01", "", ed$Array)
densityPlot(dat2, sampGroups = ed$fila)
# Es la columna 6
for (i in 1:nrow(ed)) {
  ed$columna[i] = substr(ed$Slide[i], 10, 12)
}
densityPlot(dat2, sampGroups = ed$columna)
# Es la de sufijo 053
# ES LA MUESTRA PAT_187 --> la quitamos 


# QC EWAS tools report
all_metrics = do.call("rbind", lapply(ed$Basename, function(x){
  f = read_idats(idat_files = x) 
  metrics = control_metrics(f)
  unlist(metrics)
}))
rownames(all_metrics) = ed$Sample_Name
thresholds = data.frame(variable = colnames(all_metrics), 
                        value = c(0, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5 , 5))


#Se obtiene el gráfico de los parámetros de calidad del paquete Ewastools
plots = lapply(thresholds$variable, function(x){
  toplot = data.frame(values = all_metrics[, x],
                      aux.y = rnorm(nrow(all_metrics)), 
                      sample_name = ed$Sample_Name,
                      stringsAsFactors = F)
  # toplot$color = "CTRL"
  toplot$color[grep("PAT", toplot$sample_name)] = "PAT"
  toplot$color[grep("m", toplot$sample_name)] = "CTRL_m"
  toplot$color[grep("u", toplot$sample_name)] = "CTRL_u"
  toplot$sample_name[toplot$color == "PAT"] = ""
  ggplot(toplot, aes(x=values, y=aux.y, color = color)) +
    geom_point()+
    geom_vline(xintercept = thresholds$value[thresholds$variable==x], linetype = 2)+
    xlab("")+
    ylab("")+
    ggtitle(x)+ 
    scale_color_manual(values = c("PAT" = "#8B87BF", "CTRL_m" = "#74C4AC", "CTRL_u" = "#E2833C")) +
    theme_light()+
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=10),
          axis.text.y = element_blank(),
          plot.title = element_text(size=10, hjust = 0.5), 
          legend.position = "none")
  
})
do.call("grid.arrange", c(plots, ncol = 3))

# En el RUN 3 hay una muestra que se va en BCII --> PAT_132 la quitamos 
toplot = data.frame(values = all_metrics[, "Bisulfite Conversion II (Bkg)"],
                    aux.y = rnorm(nrow(all_metrics)), 
                    sample_name = ed$Sample_Name,
                    stringsAsFactors = F)
toplot$color = "CTRL"
toplot$color[grep("PAT", toplot$sample_name)] = "PAT"
ggplot(toplot, aes(x=values, y=aux.y, color = color, label = sample_name))+
  geom_point()+
  geom_vline(xintercept = thresholds$value[thresholds$variable=="Bisulfite Conversion II (Bkg)"], linetype = 2)+
  ylab("")+
  ggtitle("Bisulfite Conversion II (Bkg)")+
  theme_light()+
  theme(axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.text.y = element_blank(),
        plot.title = element_text(size=15, hjust = 0.5), 
        legend.position = "bottom", 
        legend.text = element_text(size = 10),
        legend.title = element_blank()
        )

# Cálculo valores Beta ---------------------------------------------------------
# Se obtienen los valores Beta a partir de los iDATs
datBetas = getBeta(dat2)

#### Media de los cpgs duplicados ####
cpgs = unlist(lapply(rownames(datBetas), function(x){
  unlist(strsplit(x, split = "_"))[1]
}))
index2rm = which(duplicated(cpgs))

dups = unique(cpgs[index2rm])
length(dups)
# 5222

# Cantidad total de sitios cpgs sin duplicados
tokeep = cpgs[!cpgs%in%dups]
length(tokeep)
# 925374

# Betas sin Duplicados
datBetas_sinDup = datBetas[cpgs%in%tokeep, ]
rownames(datBetas_sinDup) = unlist(lapply(rownames(datBetas_sinDup), function(x){
  unlist(strsplit(x, split = "_"))[1]
}))

# Calcular la media de los duplicados 
n <<- 0
betas2add = do.call("rbind", lapply(dups, function(x){
  n <<- n+1
  if (n%%100 == 0){
    cat(n, "\n")
  }
  cpg_dup = grep(x, rownames(datBetas), value = T)
  apply(datBetas[cpg_dup, ], 2, mean)
}))
rownames(betas2add) = dups

# Unir los datos de las sondas no duplicadas y la media de las duplicadas
datBetas_sinDup = rbind(datBetas_sinDup, betas2add)

# Cambiar los nombres de columna a PAT_XXX
colnames(datBetas_sinDup) = ed[colnames(datBetas_sinDup), "Sample_Name"]


# PATIENT CLINICAL INFO --------------------------------------------------------
#Cargamos el excel con los datos clinicos de cada RUN

# RUN 1 
ed_pat1 = read.delim(paste0(root, "data/RUNS/RUN_1_METH_DATOS_CLINICOS.csv"),
                       sep = ",", as.is = T, row.names = 1) 
# RUN 2 
ed_pat2 = read.delim(paste0(root, "data/RUNS/RUN_2_METH_DATOS_CLINICOS.csv"),
                       sep = ",", as.is = T, row.names = 1) 
# RUN 3
ed_pat3 = read.delim(paste0(root, "data/RUNS/RUN_3_METH_DATOS_CLINICOS.csv"),
                       sep = ",", as.is = T, row.names = 1) 
# RUN 4 
ed_pat4 = read.delim(paste0(root, "data/RUNS/RUN_4_METH_DATOS_CLINICOS.csv"),
                       sep = ",", as.is = T, row.names = 1) 

# Juntar la información 
ed_pat = rbind(ed_pat1, ed_pat2, ed_pat3, ed_pat4)
str(ed_pat)

# Unificamos la información 
ed_pat$BIRTHDAY= as.Date(ed_pat$BIRTHDAY, format = "%d/%m/%Y")
ed_pat$BIOPSY_DATE = as.Date(ed_pat$BIOPSY_DATE, format = "%d/%m/%Y")
ed_pat$BIOPSY_TIME = as.POSIXct(paste(ed_pat$BIOPSY_DATE, ed_pat$BIOPSY_TIME, sep = " "), 
                                  format = "%Y-%m-%d %H:%M")
ed_pat$TR_DATE =trimws(ed_pat$TR_DATE)
ed_pat$TR_DATE[ed_pat$TR_DATE %in% c("NA", "na", "")] <- NA 
ed_pat$TR_DATE <- as.Date(ed_pat$TR_DATE, format = "%md/%m/%Y")
ed_pat$TR_No_EMBRYOS =trimws(ed_pat$TR_No_EMBRYOS)
ed_pat$TR_No_EMBRYOS[ed_pat$TR_No_EMBRYOS %in% c("NA", "na", "")] <- NA 
ed_pat$P4_DATE <- as.Date(ed_pat$P4_DATE, format = "%d/%m/%Y")
ed_pat$P4_TIME <- as.POSIXct(paste(ed_pat$P4_DATE, ed_pat$P4_TIME, sep = " "), 
                               format = "%Y-%m-%d %H:%M")
str(ed_pat)

# Añadimos la edad continua: "AGE2"
ed_pat$AGE2 = as.numeric(difftime(ed_pat$BIOPSY_DATE, 
                                  ed_pat$BIRTHDAY, 
                                  units = "days"))/365.25

# Añadimos las horas de progesterona 
ed_pat$P4_HOURS = difftime(ed_pat$BIOPSY_TIME, ed_pat$P4_TIME, units = "hours")

# Mantenemos el orden de la matriz de betas
ed_pat = ed_pat[colnames(datBetas_sinDup), ]


# PATIENT EXPERIMENTAL INFO ----------------------------------------------------
# Cargamos el excel con los datos experimentales de la extracción de DNA de cada RUN
# RUN 1
ed_exp1 =  read.delim(paste0(root, "data/RUNS/RUN_1_METH_EXPERIMENTAL.csv"), 
                      sep = ",", as.is = T, row.names = 1)
# RUN 2
ed_exp2 =  read.delim(paste0(root, "data/RUNS/RUN_2_METH_EXPERIMENTAL.csv"), 
                      sep = ",", as.is = T, row.names = 1)
# RUN 3
ed_exp3 =  read.delim(paste0(root, "data/RUNS/RUN_3_METH_EXPERIMENTAL.csv"), 
                      sep = ",", as.is = T, row.names = 1)
# RUN 4
ed_exp4 =  read.delim(paste0(root, "data/RUNS/RUN_4_METH_EXPERIMENTAL.csv"), 
                      sep = ",", as.is = T, row.names = 1)

# Juntamos la información 
ed_exp = rbind(ed_exp1, ed_exp2, ed_exp3, ed_exp4)
str(ed_exp)

# Mantenemos el orden de la matriz de betas
ed_exp = ed_exp[colnames(datBetas_sinDup), ]

# Controles --------------------------------------------------------------------

colnames(datBetas_sinDup)[colnames(datBetas_sinDup) == "CTRL_m"] = "CTRL_m_R1"
colnames(datBetas_sinDup)[colnames(datBetas_sinDup) == "CTRL_u"] = "CTRL_u_R1"

rownames(expInformation)[rownames(expInformation) == "CTRL_m"] = "CTRL_m_R4"
rownames(expInformation)[rownames(expInformation) == "CTRL_u"] = "CTRL_u_R4"

rownames(patients)[rownames(patients) == "CTRL_m"] = "CTRL_m_R1"
rownames(patients)[rownames(patients) == "CTRL_u"] = "CTRL_u_R1"


# Quitamos las muestras que no pasan los controles de calidad ------------------
# RUN3 -> PAT_132
Betas = datBetas_sinDup[, colnames(datBetas_sinDup) != "PAT_132"]
ed_pat = ed_pat[rownames(ed_pat) != "PAT_132", ]
ed_exp = ed_exp[rownames(ed_exp) != "PAT_132", ]

# RUN4 -> PAT_187 
Betas = Betas[, colnames(Betas) != "PAT_187"]
ed_pat = ed_pat[rownames(ed_pat) != "PAT_187", ]
ed_exp = ed_exp[rownames(ed_exp) != "PAT_187", ]

# Save data --------------------------------------------------------------------

write.csv(x = Betas, file = paste0(results_folder,"matriz_betas.csv"),
          row.names = T, quote = F) 

save(Betas, ed_pat, ed_exp, 
     file = paste0(results_folder, "meth_info.RData"))






# DESCRIPCIÓN: Análisis exploratorio de los datos con limma
# Se corrige el efecto de Batch del proyecto 
# DEA

# Settings ----------------------------------------------------------------

rm(list=ls())

# Working directory
root = "/data/local/mvaquerizo/transcriptomica/"
setwd(root)

# Librerías
library(gplots)
library(pvca)
library(stringr)
library(org.Hs.eg.db)
library(NOISeq)

# Funciones 
source("Funciones_R.R")

# Carpeta de resultados
results_folder = paste0(root, "results/06_exploratory_and_DEA_limma/")
dir.create(results_folder)

# Carga de datos ---------------------------------------------------------------

load(paste0(root, "results/05_select_transcriptome_samples/conteos_finales_E_R.RData"), verbose = T)
# Loading objects:
#   pat_final_counts_filt
#   metadata_pat

str(pat_final_counts)
str(metadata_pat)

# Generamos paleta de color envejecida - rejuvenecida 
palete_g_acc = c("#3FBCC3", "#B75180")

# Voom + quantiles (limma) ------------------------------------------------

#### Filtro bajo conteos ####
# Filtra por bajos conteos, pero respetando diferencias entre grupos

zero_genes = rownames(pat_final_counts)[apply(pat_final_counts,1,sum) == 0]
to_keep = rownames(pat_final_counts)[apply(pat_final_counts,1,sum) != 0]
pat_final_counts = pat_final_counts[to_keep, ]
counts_filtered = filtered.data(dataset = pat_final_counts,
                                factor = metadata_pat$g_acc,
                                norm = F,
                                method = 1,
                                cpm = 1.5)
# Nos quedamos con 13094 genes


#### Normalización ####
voom_obj = voom(counts = counts_filtered, normalize.method = "quantile", plot = T)
counts_norm = voom_obj$E

# Boxplot sin normalizar
plot_Boxplot(counts_filtered, ed = metadata_pat, condition = "g_acc", 
             title = "Datos sin normalizar", colors = palete_g_acc)
# Boxplot normalizado
plot_Boxplot(counts_norm, ed = metadata_pat, condition = "g_acc", 
             title = "Datos normalizados", colors = palete_g_acc) 


#### Batch effect ####
summary(metadata_pat$DV200)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 54.96   81.39   84.12   83.27   85.81   90.45 
metadata_pat$g_DV200 = cut(metadata_pat$DV200, breaks = c(54, 81.7, 84.4, 85.8, 91))
table(metadata_pat$g_DV200, useNA = "ifany")
# (54,81.4] (81.4,84.2] (84.2,85.8]   (85.8,91] 
#        18          18          20          16

summary(metadata_pat$RIN)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   1.000   4.700   6.250   6.121   7.625   9.100       1
metadata_pat$g_RIN = cut(metadata_pat$RIN, breaks = c(0, 4.6, 6.2, 7.8, 10))
table(metadata_pat$g_RIN, useNA = "ifany")
# (0,4.7] (4.7,6.2] (6.2,7.6]  (7.6,10]      <NA> 
#      18        19        16        18         1 

pData = AnnotatedDataFrame(metadata_pat[, c("g_DV200", "g_RIN", "project")])
datset = ExpressionSet(counts_norm, phenoData = pData)

# PVCA
pvca_res = pvcaBatchAssess(abatch = datset,
                           batch.factors = c("g_DV200", "g_RIN", "project"),
                           threshold = 1)

to_plot <-  data.frame(WC = t(pvca_res$dat), feature = pvca_res$label)
to_plot$feature <- str_replace_all(to_plot$feature, c("g_RIN" = "RIN", 
                                                      "g_DV200" = "DV200", 
                                                      "project" = 'Proyecto',
                                                      "resid"= "Residual"))
to_plot$feature = factor(to_plot$feature, levels = to_plot$feature[order(to_plot$WC)])
to_plot$color <- "inter"
to_plot[grep(":",to_plot$feature, fixed = T, invert = T),"color"] = "variable"
to_plot[grep("Residual",to_plot$feature),"color"] = "residual"

ggplot(data = to_plot[to_plot$color!="inter",], aes(x=feature, y = WC*100, fill = color)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(WC * 100, 2), "%")), vjust = -0.2, size = 5) +
  xlab("")+
  ylab ("Variabilidad (%)") +
  ggtitle("Análisis PVCA (pre)") +
  geom_hline(yintercept = 10, linetype = 2, color = "darkgrey")+
  scale_fill_manual(values = c(inter = "darkgoldenrod3",
                               variable = "blue3",
                               resid = "black")) +
  default_theme() + 
  theme(x.axis.angle = 45)

ggsave(paste0(results_folder, "pvca_pre.jpg"), dpi = 500)

# PCA
plot_PCAscores(dat = counts_norm, ed = metadata_pat, condition1 = "project", 
               title = "PCA (proyecto)", colors = c("#66CDAA", "#CD9ABC")) + 
  coord_equal()


# Quitamos efecto de batch de proyecto
counts_norm_remove = limma::removeBatchEffect(x=counts_norm, 
                                              batch = metadata_pat$project)

# PVCA con bath corregido
pData2 = AnnotatedDataFrame(metadata_pat[, c("g_DV200", "g_RIN", "project")])
datset2 = ExpressionSet(counts_norm_remove, phenoData = pData2)
pvca_res2 = pvcaBatchAssess(abatch = datset2,
                           batch.factors = c("g_DV200", "g_RIN", "project"),
                           threshold = 1)

to_plot =  data.frame(WC = t(pvca_res2$dat), feature = pvca_res2$label)

to_plot$feature = str_replace_all(to_plot$feature, c("g_RIN" = "RIN",
                                                     "g_DV200" = "DV200", 
                                                     "project" = "Proyecto",
                                                     "resid"= "Residual"))
to_plot$feature = factor(to_plot$feature, levels = to_plot$feature[order(to_plot$WC)])
to_plot$color = "inter"
to_plot[grep(":",to_plot$feature, fixed = T, invert = T),"color"] = "variable"
to_plot[grep("Residual",to_plot$feature),"color"] = "residual"


ggplot(data = to_plot[to_plot$color!="inter",], aes(x=feature, y = WC*100, fill = color)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(WC * 100, 2), "%")), vjust = -0.2, size = 5) +
  xlab("")+
  ylab ("Variabilidad (%)") +
  ggtitle("Análisis PVCA (post)") +
  geom_hline(yintercept = 10, linetype = 2, color = "darkgrey")+
  scale_fill_manual(values = c(inter = "darkgoldenrod3",
                               variable = "blue3",
                               resid = "black")) +
  default_theme() + 
  theme(x.axis.angle = 45)

ggsave(paste0(results_folder, "pvca_post.jpg"), dpi = 500)

plot_PCAscores(dat = counts_norm_remove, ed = metadata_pat2, condition1 = "project", 
               title = "PCA proyecto (batch corregido)", colors = c("#66CDAA", "#CD9ABC"))  + 
  coord_equal()


# DEA ---------------------------------------------------------------------

dea_res = diffExprAnalysis(dat = counts_norm_remove, ed = metadata_pat, condition = "g_acc")
head(dea_res$`E-R`)

sigs = dea_res$`E-R`[dea_res$`E-R`$adj.P.Val<=0.05, ]
dim(sigs)
# 0  7


# Save ---------------------------------------------------------------------

save(counts_norm_remove, dea_res, file = paste0(results_folder, "DEA_E_R_limma.RData"))
write.csv(counts_norm_remove, file = paste0(results_folder, "conteos_finales.csv"))
write.csv(metadata_pat, file = paste0(results_folder, "metadata.csv"))

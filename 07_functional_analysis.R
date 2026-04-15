# mar.vaquerizo@ivirma.com 

# DESCRIPTION: Análisis funcional 
# GSEA con clusterprofiler ordenando por logFC
# Anotación con GO y KEGG

# Settings ----------------------------------------------------------------

rm(list=ls())

# Working directory
root = "/data/local/mvaquerizo/transcriptomica/"
setwd(root)


# Library
library(gplots)
library(stringr)
library(tibble)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
library(KEGG.db)

# Funciones 
source("../TFM/Funciones.R")

# Carpeta de resultados
results_folder = paste0(root, "results/07_functional_analysis/")
dir.create(results_folder)

# Carga de datos ---------------------------------------------------------------

load(paste0(root, "results/05_select_transcriptome_samples/conteos_finales_E_R.RData"), verbose = T)
# Loading objects:
#   pat_final_counts_filt
#   metadata_filt

load(paste0(root, "results/06_exploratory_and_DEA_limma/DEA_E_R.RData"), verbose = T)
# Loading objects:
#   final_counts
#   pat_final_counts
#   metadata


# GSEA -------------------------------------------------------------------------

#### GO ####
res = dea_res$`E-R`
# Sacamos los genes y ordenamos por logFC decreciente
geneList = res$logFC 
names(geneList) = rownames(res)  
geneList = sort(geneList, decreasing = T)  
head(geneList)

gsea_GO = gseGO(geneList,
                 ont = "ALL",
                 OrgDb = org.Hs.eg.db, 
                 keyType = "SYMBOL", 
                 minGSSize = 10, 
                 maxGSSize = 300,
                 seed = 123)
res_gsea_go = gsea_GO@result
write.csv(res_gsea_go, file = paste0(results_folder, "funciones_gsea_GO_E_R.csv"))


#### KEGG ####

# Sacamos los entrez id 
entrez_ids = mapIds(org.Hs.eg.db, keys=names(geneList), column="ENTREZID", keytype="SYMBOL")
df = data.frame("symbol" = names(entrez_ids), "entrez" = entrez_ids)
all(df$symbol == names(geneList))
# TRUE 

# Ordenamos los genes (entrez) por logFC decreciente
df$logFC = geneList
df_filt = df %>% 
  arrange(desc(logFC)) %>% 
  filter(!(is.na(entrez))) %>% 
  distinct(entrez, .keep_all = T)

geneList_kegg = df_filt$logFC
geneList_kegg = sort(geneList_kegg, decreasing = T)
names(geneList_kegg) = as.character(df_filt$entrez)

gsea_KEGG = gseKEGG(geneList_kegg, 
                    organism = "hsa", 
                    keyType = 'ncbi-geneid',
                    minGSSize = 10, 
                    maxGSSize = 300, 
                    use_internal_data = T,
                    seed = 123)
res_gsea_kegg = gsea_KEGG@result
write.csv(res_gsea_kegg, file = paste0(results_folder, "funciones_gsea_KEGG_E_R.csv"))

#### Unir los resultados ####
# Combinar los resultados de GO y KEGG
all(colnames(res_gsea_go) == colnames(res_gsea_kegg))
# FALSE 
# En los resultados de KEGG no está la columna ontology, la añadimos
res_gsea_kegg$ONTOLOGY = "KEGG"
res_gsea_kegg = res_gsea_kegg[, colnames(res_gsea_go)]
all(colnames(res_gsea_go) == colnames(res_gsea_kegg))
# TRUE
res_gsea_bind = rbind(res_gsea_go, res_gsea_kegg)
write.csv(res_gsea_bind, file = paste0(results_folder, "funciones_gsea_bind_E_R.csv"))


#### Merge con atlas funciones ####
atlas = read.csv("data/proyectos/functional_groups_atlas.csv", row.names = F)
res_gsea_bind$Description = tolower(res_gsea_bind$Description)
funciones = merge(res_gsea_bind[, 1:5], atlas[, 2:ncol(atlas)], by = "Description",
                  all.x = T, all.y = F)
write.csv(funciones, file = paste0(results_folder, "funciones_gsea_atlas.csv"))



# Interpretación grupos funcionales --------------------------------------------
res_gsea_bind = read.csv("results/04_functional_analysis/funciones_gsea_bind_E_R.csv")
grupos_funciones = read.csv("results/04_functional_analysis/funciones_E_R_agrupadas_final.csv")
all(grupos_funciones$Description %in% res_gsea_bind$Description)
# TRUE 
grupos_funciones_gsea = merge(grupos_funciones[, c("Description", "Grupos_finales")],
                              res_gsea_bind, 
                              by = "Description")
grupos_funciones_gsea$updown = ifelse(grupos_funciones_gsea$NES < 0, "Down-regulated", "Up-regulated")

table(grupos_funciones_gsea$Grupos_finales, grupos_funciones_gsea$updown)
#                                               Down-regulated Up-regulated
# Actividad ribosomal                                       14            0
# Adhesion                                                   4           18
# Desarrollo de tejidos                                      0           13
# Deteccion de estimulos y percepcion sensorial             13            0
# Epigenetica y regulacion genica                            0            7
# Estres oxidativo y detoxificacion                         17            0
# Hormonas y señalizacion                                    1           14
# Inflamacion e inmunidad innata                            64            0
# Inmunidad adaptativa                                      42            0
# Metabolismo glucidos y energia                             0           27
# Metabolismo lipidos                                        0           20
# Otros                                                      4           10
# Proteostasis y glicosilacion                               0           21
# Regulacion de iones                                        0           34
# Regulacion respuesta inmune                               44            0
# Sistema nervioso                                           0           18
# Transporte celular                                         0           22

write.csv(grupos_funciones_gsea, file = paste0(results_folder, "funciones_E_R.csv"))


# Gráfico de pie de las proporciones de cada grupo de funciones 
# Generamos df con los grupos, el número de funciones up-down y el color
# rojo cuando todas up, azul todas down y gris alguna up-down
tabla = table(grupos_funciones_gsea$Grupos_finales, grupos_funciones_gsea$updown)
toplot = as.data.frame.matrix(tabla)
toplot$Grupo = rownames(toplot)
toplot = toplot %>% 
  mutate(N = `Down-regulated` + `Up-regulated`)
toplot$color = ifelse(toplot$`Down-regulated` == 0 & toplot$`Up-regulated` > 0, "#FF8288", 
                      ifelse(toplot$`Down-regulated` > 0 & toplot$`Up-regulated` == 0, "#B2D0E2", "gray70"))
toplot["Sistema inmune", "color"] = "#B2D0E2" # como solo hay 1 up lo pongo azul
# Ordenamos para que salgan juntos por colores
toplot = toplot %>% 
  arrange(toplot, color) %>% 
  mutate(Grupo = factor(Grupo, levels = Grupo))
colores = setNames(toplot$color, toplot$Grupo)
# pie-chart
ggplot(toplot, aes(x=1, y = N, fill = Grupo)) + 
  geom_bar(stat = "identity", width = 1, color = "white") +
  geom_text(aes(x = 1.3, label = N), position = position_stack(vjust = 0.5), 
            size = 6) +
  coord_polar(theta = "y", start = 0) + 
  theme_minimal() +
  theme_void() + 
  scale_fill_manual(values = colores)

ggsave("plots/pie_char_funciones.jpg", dpi = 500, width = 11, height = 7)

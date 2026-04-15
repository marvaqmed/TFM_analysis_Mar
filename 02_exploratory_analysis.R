# DESCRIPCIÓN: Análisis exploratorio de los datos
# Se eliminan los NAs de la matriz de betas 
# Se comprueba el agrupamiento de los controles de metilación 
# Detección de outliers
# Comparar características clínicas y demográficas de la población control y FIV

# Settings ----------------------------------------------------------------

rm(list=ls())

# Working directory
root = "/data/local/mvaquerizo/metilacion/"
setwd(root)

# Librerías
library(gridExtra)
library(ggplot2)
library(pvca)
library(Biobase)
library(ggforce)
library(HotellingEllipse)
library(FactoMineR)
library(tibble)

# Funciones 
source("../TFM/Funciones.R")

# Carpeta de resultados
results_folder = paste0(root, "results/02_exploratory_analysis/")
dir.create(results_folder)


# Carga de datos ---------------------------------------------------------------
load("results/01_fromRAW_to_beta/meth_info.RData", verbose = T)
# Loading objects:
# Betas
# ed_pat
# ed_exp

# Eliminar los NAs de la matriz de Betas ---------------------------------------
# Encontrar los Betas de las muestras que tengan "NA"
apply(Betas, 2, function(x){sum(is.na(x))})
dim(Betas)
# 930596    190

#Se eliminan los Betas que tengan "NA"
index2remove = apply(Betas, 1, function(x){any(is.na(x))})
sum(index2remove)
# 2286

datBetas2pcas = Betas[!index2remove, ]
dim(datBetas2pcas)
# 928310    190

# Análisis eploratorio ---------------------------------------------------------

#### Comportamiento control metilado positivo y negativo ####
ed_exp$control = "muestra"
ed_exp[grep("CTRL_m", rownames(ed_exp)), "control"] = "control_m"
ed_exp[grep("CTRL_u", rownames(ed_exp)), "control"] = "control_u"
ed_exp$control = factor(ed_exp$control, levels = c("muestra", "control_m", "control_u"))
ed_exp$muestra = rownames(ed_exp)
plot_PCAscores(dat = datBetas2pcas, ed = ed_exp, condition1 = "control") +
  coord_equal() + 
  ggtitle("Controles metilados y no metilados") +
  scale_color_manual(values = c("#8B87BF", "#74C4AC", "#E2833C")) +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.text = element_text(size = 15), 
        plot.title = element_text(size=15, hjust = 0.5)) 


# Eliminamos los controles 
all.equal(colnames(datBetas2pcas), rownames(ed_exp))
# TRUE
datBetas2pcas = datBetas2pcas[, ed_exp$control == "muestra"]
Betas = Betas[, ed_exp$control == "muestra"]
ed_pat = ed_pat[ed_exp$control == "muestra", ]
ed_exp = ed_exp[ed_exp$control == "muestra", ]


#### Outliers ####
# Calculamos el PCA y nos quedamos con las 2 primeras componentes
# Uso de elipses de confianza basadas en distancias de Hotelling
pca = PCA(t(datBetas2pcas))
pca_scores = as_tibble(pca$ind$coord)
res_2PCs <- ellipseParam(pca_scores, k = 2, pcx = 1, pcy = 2)
str(res_2PCs)
a1 = res_2PCs$Ellipse$a.95pct
b1 = res_2PCs$Ellipse$b.95pct
T2 = res_2PCs$Tsquare$value

coord_2PCs_95 <- ellipseCoord(pca_scores, pcx = 1, pcy = 2, conf.limit = 0.99, pts = 500)

# Con labels
ggplot(pca_scores, aes(x = Dim.1, y = Dim.2)) +
  geom_path(data = coord_2PCs_95, aes(x, y), color = "cornsilk4") +
  geom_point(aes(fill = ed_exp$intrarun, color=ed_exp$intrarun), shape = 21, size = 3) +
  geom_text_repel(aes(label = ed_exp$muestra, color = ed_exp$intrarun),
                  size = 4, point.padding = unit(0.2, "lines"), segment.size = 0)+
  labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC2", 
       x = paste0("PC1:", round(pca$eig[1,2], 1), "%"), 
       y = paste0("PC2:", round(pca$eig[2,2], 1), "%")) +
  scale_color_manual("#8B87BF") + 
  default_theme() + 
  coord_equal() + 
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=18),
        legend.position = "none", 
        plot.title = element_text(size=20, hjust = 0.5))

# Sin labels
ggplot(pca_scores, aes(x = Dim.1, y = Dim.2)) +
  geom_path(data = coord_2PCs_95, aes(x, y), color = "cornsilk4") +
  geom_point(aes(fill = ed_exp$control, color=ed_exp$control), shape = 21, size = 3) +
  labs(title = "Detección de outliers", 
       x = paste0("PC1:", round(pca$eig[1,2], 1), "%"), 
       y = paste0("PC2:", round(pca$eig[2,2], 1), "%")) +
  scale_fill_manual(values = "#8B87BF") + 
  scale_color_manual(values = "#8B87BF") +
  default_theme() + 
  coord_equal() +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=18),
        legend.position = "none", 
        plot.title = element_text(size=20, hjust = 0.5))
  

# Eliminamos outliers
outliers = c("PAT_41")
ed_exp = ed_exp[!(rownames(ed_exp) %in% outliers), ]
dim(ed_exp)
# 169  21
ed_pat = ed_pat[!(rownames(ed_pat) %in% outliers), ]
dim(ed_pat)
# 169   21
Betas = Betas[, !(colnames(Betas) %in% outliers)]
dim(Betas)
# 930596    169
datBetas2pcas = datBetas2pcas[, !(colnames(datBetas2pcas) %in% outliers)]
dim(datBetas2pcas)
# 928310    169



# Comportamiento experimental variables ----------------------------------------

#### Kit DNA #####
table(ed_exp$DNA_extraction._kit, useNA = "ifany")
# DNA_QIAamp   DNA/RNA/miRNA 
#        140              29 

plot_PCAscores(dat = datBetas2pcas, ed = ed_exp, condition1 = "DNA_extraction._kit") +
  coord_equal() +
  ggtitle("Kit de extracción ADN") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))
  

#### A260/280 #####
summary(ed_exp$DNA_A260.280)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1.750   1.870   1.880   1.885   1.890   2.060       2 
# Creamos nueva columna agrupando por intervalos
ed_exp$DNA_A260.280_gr = ifelse(ed_exp$DNA_A260.280 < 1.87, "<1.87", 
                                ifelse(ed_exp$DNA_A260.280 < 1.88, "1.87-1.88", 
                                       ifelse(ed_exp$DNA_A260.280 < 1.89, "1.88-1.89", ">1.89")))

plot_PCAscores(dat = datBetas2pcas, ed = ed_exp, condition1 = "DNA_A260.280_gr") +
  coord_equal() +
  ggtitle("A260/A280") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))


#### A260/230 #####
summary(ed_exp$DNA_A260.230)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.340   1.270   1.670   1.612   2.005   2.490       2 
# Creamos nueva columna agrupando por intervalos
ed_exp$DNA_A260.230_gr = ifelse(ed_exp$DNA_A260.230 < 1.3, "<1.3", 
                                ifelse(ed_exp$DNA_A260.230 < 1.7, "1.3-1.7", 
                                       ifelse(ed_exp$DNA_A260.230 < 2, "1.7-2.0", ">2.0")))
ed_exp$DNA_A260.230_gr = factor(ed_exp$DNA_A260.230_gr, 
                                levels = c("<1.3", "1.3-1.7", "1.7-2.0", ">2.0"))
table(ed_exp$DNA_A260.230_gr, useNA = "ifany")
# <1.3 1.3-1.7 1.7-2.0    >2.0    NA's
#   45      43      38      41       2 
plot_PCAscores(dat = datBetas2pcas, ed = ed_exp, condition1 = "DNA_A260.230_gr") + 
  coord_equal() +
  ggtitle("A260/A230") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))

#### RUN metilación #####
table(ed_exp$RUN.METH, useNA = "ifany")
# PIL_EA_1 PIL_EA_2 PIL_EA_3 PIL_EA_4 
#       45       45       45       34
plot_PCAscores(dat = datBetas2pcas, ed = ed_exp, condition1 = "RUN.METH") + 
  coord_equal() +
  ggtitle("Run secuenciación") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))


# Comportamiento de variables clínicas -----------------------------------------

#### BMI ####
summary(ed_pat$BMI)
# Creamos nueva columna agrupando por intervalos
ed_pat$BMI_g = cut(ed_pat$BMI, breaks = c(17, 21, 22.3, 24.5, 35))
ed_pat$BMI_g = factor(ed_pat$BMI_g, levels = c("(17,21]", "(21,22.3]", "(22.3,24.5]", "(24.5,35]"),
                      labels = c("<21", "21-22.3", "22.3-24.5", ">24.5"))
table(ed_pat$BMI_g)
# 17 - 21  21 - 22,3, 22,3 - 24,5   24,5 - 35 
#      46          39          41          40 

plot_PCAscores(dat = datBetas2pcas, ed = ed_pat, condition1 = "BMI_g") + 
  coord_equal() + 
  ggtitle("BMI") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))

#### Oocyte origin ####
ed_pat$TR_OOCYTE_ORIGIN = factor(ed_pat$TR_OOCYTE_ORIGIN, levels = c("Donated", "Own", NA), 
                                 labels = c("Donado", "Propio"))
plot_PCAscores(dat = datBetas2pcas, ed = ed_pat, condition1 = "TR_OOCYTE_ORIGIN") +
  coord_equal() + 
  ggtitle("Origen del ovocito") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))

#### Horas de P4 ####
summary(as.numeric(ed_pat$P4_HOURS))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 82.0   111.5   116.2   117.3   122.9   148.6      11 
# Creamos nueva columna agrupando por intervalos
ed_pat$P4_hours_g = cut(as.numeric(ed_pat$P4_HOURS), breaks = c(81, 111.5, 116.5, 123, 150))
table(ed_pat$P4_hours_g)
# (81,112] (112,116] (116,123] (123,150] 
# 41        40        38        39
ed_pat$P4_hours_g = factor(ed_pat$P4_hours_g, levels = c("(81,112]","(112,116]","(116,123]","(123,150]"),
                           labels = c("<112", "112-116", "116-123", ">123"))
plot_PCAscores(dat = datBetas2pcas, ed = ed_pat, condition1 = "P4_hours_g") +
  coord_equal() + 
  ggtitle("Horas de progesterona") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))

#### Patología uterina benigna #### 
table(ed_pat$BUP)
# FALSE  TRUE 
#   101    45 
ed_pat$BUP = ifelse(ed_pat$BUP == TRUE, "BUP", "No BUP")
plot_PCAscores(dat = datBetas2pcas, ed = ed_pat, condition1 = "BUP") +
  coord_equal() + 
  ggtitle("Presencia de patologías uterinas benignas (BUP)") + 
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        legend.position = "bottom", 
        plot.title = element_text(size=20, hjust = 0.5))


# PVCA --------------------------------------------------------------------

# Añadimos como variable de interés el proyecto del grupo del que proceden las muestras
ed_exp$Proyecto = ifelse(ed_exp$STUDY == "REPROCYCLE", "REPROCYCLE", 
                        ifelse(ed_exp$STUDY %in% c('TED', 'PIF'), 'PIF', 'ENCLIVA'))

pData =  AnnotatedDataFrame(cbind(ed_exp[, c("DNA_extraction._kit", "DNA_A260.230", 
                                             "DNA_A260.280", "RUN.METH", "Proyecto")]))
dat_set = ExpressionSet(datBetas2pcas, phenoData = pData, byrow=FALSE)
pvca_res = pvcaBatchAssess(abatch = dat_set, 
                           batch.factors = c("DNA_extraction._kit", "DNA_A260.230", "DNA_A260.280", "RUN.METH", "Proyecto"), 
                           threshold = 1)

# Representación gráfica
toplot = data.frame(WC = t(pvca_res$dat) * 100,
                    class = pvca_res$label, 
                    stringsAsFactors = F)
toplot$class = gsub("RUN.METH", "Tanda", toplot$class)
toplot$class = gsub("DNA_extraction._kit", "Kit_DNA", toplot$class)
toplot$class = gsub("DNA_A260.230", "A260_A230", toplot$class)
toplot$class = gsub("DNA_A260.280", "A260_A280", toplot$class)
toplot$class[toplot$class == "resid"] = "Residual"
toplot$class = factor(toplot$class, levels = toplot$class[order(toplot$WC)])
toplot$color = "inter"
toplot$color[grep(":", toplot$class, fixed = T, invert = T)] = "variable"
toplot$color[toplot$class=="Residual"] = "residual"
ggplot(data = toplot[toplot$color!="inter",], aes(x=class, y = WC, fill = color)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(WC, 2), "%")), vjust = -0.2, size = 5) +
  xlab("")+
  ylab ("Variabilidad (%)") +
  geom_hline(yintercept = 10, linetype = 2, color = "darkgrey")+
  scale_fill_manual(values = c(inter = "darkgoldenrod3",
                               variable = "blue3",
                               resid = "black")) +
  default_theme() + 
  theme(x.axis.angle = 45) +
  ggtitle('Analisis PVCA')


# Tabla demográfica  ------------------------------------------------------

# Sacamos la tabla demográfica para comparar las variables clínicas y demográficas
# de las mujeres del grupo de pacientes y de control
demo_dat = ed_pat[, c("AGE2", "BMI")]
demo_dat$P4_HOURS = as.numeric(ed_pat$P4_HOURS)
demo_dat$grupo = ifelse(ed_pat$STUDY == "REPROCYCLE", "Control", "Paciente")
tabla = demographic_table(demo_dat, demo_dat$grupo)


# Save --------------------------------------------------------------------

save(Betas, ed_exp, ed_pat, file = paste0(results_folder, "data_filtered.RData"))
save(pvca_res, file = paste0(results_folder, "PVCA_results_allcpgs.RData"))



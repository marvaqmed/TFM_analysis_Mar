# DESCRIPCIÓN: 
# Predecir la edad biológica del grupo control y de pacientes FIV
# Representación de edad biológica frente a edad cronológica 
# Cálculo aceleraciones

# Settings ----------------------------------------------------------------

rm(list=ls())

# Working directory
root = "/data/local/mvaquerizo/metilacion/"
setwd(root)

# Librerías
library(gridExtra)
library(ggplot2)
library(Biobase)
library(methylclock) 
library(VGAM)

# Funciones 
source("../TFM/Funciones.R")

# Capeta de resultados
results_folder = paste0(root, "results/03_biological_age_analysis/")
dir.create(results_folder)


# Carga de datos ---------------------------------------------------------------

load("results/02_exploratory_analysis/data_filtered.RData", verbose = T)
# Loading objects:
#   Betas
#   ed_exp
#   ed_pat


# Separar información control y FIV ---------------------------------------

control = rownames(ed_pat[ed_pat$STUDY == 'REPROCYCLE', ])
length(control)
# 19

# Seleccionamos solo la info de controles
betas_c = Betas[, control]
ed_pat_c = ed_pat[control, ]
ed_exp_c = ed_exp[control,]

# Seleccionamos la info de las pacientes (!control)
betas_pat = Betas[, !colnames(Betas)%in%control]
ed_pat_pat = ed_pat[!rownames(ed_pat)%in%control, ]
ed_exp_pat = ed_exp[!rownames(ed_exp)%in%control, ]


# Aplicación relojes epigenéticos ----------------------------------------------

clocks = c("Horvath", "BLUP", "EN")

# Predecir edad biológica controles
bioage_c = DNAmAge(x = betas_c, clocks = clocks, toBetas = F, cell.count = F, 
                 age = ed_pat_c$AGE2, normalize = T)

# Predecir edad biológica pacientes FIV 
bioage_pat = DNAmAge(x = betas_pat, clocks = clocks, toBetas = F, cell.count = F, 
                     age = ed_pat_pat$AGE2, normalize = T)

# AgeAcc: difference between DNAmAge and chronological age
# AgeAcc2: residuals obtained after regressing chronological age and DNAmAge

bioage_c = as.data.frame(bioage_c)
bioage_pat = as.data.frame(bioage_pat)
#Cambiamos los nombres de fila por los identificadores de muestras
rownames(bioage_c) = bioage_c$id
rownames(bioage_pat) = bioage_pat$id


#### AltumAge ####
# Preparar betas para Pyaging
## Controles
# Se necesita la matriz transpuesta con muestras en filas y CpG en columnas 
t_betas_c = t(betas_c)
# Eliminar NA 
col2rm = which(colSums(is.na(t_betas_c)) > 0)
col2rm = names(col2rm)
length(col2rm)
# 200
sinNA_t_betas_c = t_betas_c[, !(colnames(t_betas_c)%in%col2rm)]
dim(sinNA_t_betas_c)
# 19   930936

# Guardar info en csv 
write.csv(sinNA_t_betas_c, file = paste0(results_folder, "betas_control_sinNA.csv"))

## Pacientes 
# Se necesita la matriz transpuesta con muestras en filas y CpG en columnas 
t_betas_pat = t(betas_pat)
# Eliminar NA 
col2rm = which(colSums(is.na(t_betas_pat)) > 0)
col2rm = names(col2rm)
length(col2rm)
# 1886
sinNA_t_betas_pat = t_betas_pat[, !(colnames(t_betas_pat)%in%col2rm)]
dim(sinNA_t_betas_pat)
# 150   928715

# Guardar info en csv 
write.csv(sinNA_t_betas_pat, file = paste0(results_folder, "betas_pat_sinNA.csv"))

# Cargar info de los resultados obtenidos en Python 
pyaging_c = read.csv(paste0(root, "results/03_biological_age_analysis/bioage_control_pyaging.csv"), 
                     row.names = 1)
pyaging_pat = read.csv(paste0(root, "results/03_biological_age_analysis/bioage_pat_pyaging.csv"), 
                       row.names = 1)

#### Combinar python y R ####
all(rownames(bioage_c) == rownames(pyaging_c))
# TRUE
bioage_c = cbind(bioage_c, pyaging_c$altumage)
colnames(bioage_c)[ncol(bioage_c)] = "AltumAge"
write.csv(x = bioage_c, file = paste0(results_folder, "bioage_control_Horvath_Zhang_Altum.csv"))

all(rownames(bioage_pat) == rownames(pyaging_pat))
# TRUE
bioage_pat = cbind(bioage_pat, pyaging_pat$altumage)
colnames(bioage_pat)[ncol(bioage_pat)] = "AltumAge"
write.csv(x = bioage_pat, file = paste0(results_folder, "bioage_pat_Horvath_Zhang_Altum.csv"))


# Análisis edad biológica controles --------------------------------------------

# Seleccionamos solo las columnas de los relojes
filt_bioage_c = bioage_c[, c("Horvath", "BLUP", "EN", "AltumAge")]

#### Edad biológica vs edad cronológica #### 
plots = lapply(colnames(filt_bioage_c), function(x){
  toplot = data.frame(cronologica = ed_pat_c$AGE2, biologica = filt_bioage_c[, x])
  lims = c(min(toplot, na.rm = T), max(toplot, na.rm = T))
  ggplot(toplot, aes(x=cronologica, y=biologica))+
    geom_point(alpha = 0.6)+
    xlab("Chronological Age")+
    ylab("Biological Age")+
    ggtitle(x)+
    geom_abline(intercept = 0, linetype = 2, color = "grey25")+
    geom_smooth(method="lm", se=F)+
    default_theme()+
    theme(plot.title = element_text(size = 15, hjust=0.5),
          axis.text = element_text(size=10),
          axis.title = element_text(size=15),) + 
    coord_equal(xlim = lims, ylim = lims)
})
do.call("grid.arrange", c(plots, ncol = 4))


#### Análisis de correlación ####
corrs = do.call("rbind", lapply(colnames(filt_bioage_c), function(x){
  res = cor.test(ed_pat_c$AGE2, filt_bioage_c[, x])
  data.frame(correlation = res$estimate, 
             p_value = res$p.value)
}))
rownames(corrs) = colnames(filt_bioage_c)
corrs

# Media correlación controles 
mean(corrs$correlation)
# 0.8515271

#Aquí se dibuja la linea de correlacion en las gráficas para poder visualizarla
cor_plots = lapply(colnames(filt_bioage_c), function(x){
  toplot = data.frame(cronologica = ed_pat_c$AGE2, biologica = filt_bioage_c[, x])
  lims = c(min(toplot), max(toplot))
  ggplot(toplot, aes(x=cronologica, y=biologica))+
    geom_point(alpha = 0.6)+
    xlab("Edad cronologica")+
    ylab("Edad epigenetica")+
    ggtitle(x)+
    geom_smooth(method="lm", color = "#31a7f6", linewidth = 1.5)+
    geom_abline(intercept = 0, linetype = 2, color = "grey25")+
    default_theme()+
    theme(plot.title = element_text(size = 25, hjust=0.5),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20)) + 
    coord_equal(xlim = lims, ylim = lims)
})
do.call("grid.arrange", c(cor_plots, ncol = 4))


# Análisis edad biológica pacientes ----------------------------------------

# Seleccionamos solo las columnas de los relojes
filt_bioage_pat = bioage_pat[, c("Horvath", "BLUP", "EN", "AltumAge")]

#### Edad biológica vs edad cronológica #### 
plots = lapply(colnames(filt_bioage_pat), function(x){
  toplot = data.frame(cronologica = ed_pat_pat$AGE2, biologica = filt_bioage_pat[, x])
  lims = c(min(toplot, na.rm = T), max(toplot, na.rm = T))
  ggplot(toplot, aes(x=cronologica, y=biologica))+
    geom_point(alpha = 0.6)+
    xlab("Chronological Age")+
    ylab("Biological Age")+
    ggtitle(x)+
    geom_abline(intercept = 0, linetype = 2, color = "grey25")+
    geom_smooth(method="lm", se=F, color = "#FF7F50")+
    default_theme()+
    theme(plot.title = element_text(hjust=0.5)) + 
    coord_equal(xlim = lims, ylim = lims)
})
do.call("grid.arrange", c(plots, ncol = 4))


#### Análisis de correlación ####
corrs = do.call("rbind", lapply(colnames(filt_bioage_pat), function(x){
  res = cor.test(ed_pat_pat$AGE2, filt_bioage_pat[, x])
  data.frame(correlation = res$estimate, 
             p_value = res$p.value)
}))
rownames(corrs) = colnames(filt_bioage_pat)
corrs

# Media correlación pacientes 
mean(corrs$correlation)
# 0.5979949

#Aquí se dibuja la linea de correlacion en las gráficas para poder visualizarla
cor_plots = lapply(colnames(filt_bioage_pat), function(x){
  toplot = data.frame(cronologica = ed_pat_pat$AGE2, biologica = filt_bioage_pat[, x])
  lims = c(min(toplot), max(toplot))
  ggplot(toplot, aes(x=cronologica, y=biologica))+
    geom_point(alpha = 0.5)+
    xlab("Edad cronologica")+
    ylab("Edad biologica") +
    geom_smooth(method="lm", linewidth = 1.5, color = "#EB9050")+
    geom_abline(intercept = 0, linetype = 2, color = "grey25")+
    default_theme()+
    theme(plot.title = element_text(size = 25, hjust=0.5),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20)) + 
    coord_equal(xlim = lims, ylim = lims)
})
do.call("grid.arrange", c(cor_plots, ncol = 4))


# Pactientes FIV frente a controles --------------------------------------------

toplot <- data.frame(
  cronologica = c(ed_pat_c$AGE2, ed_pat_pat$AGE2),
  biologica = c(filt_bioage_c$EN, filt_bioage_pat$EN),
  pacientes = factor(c(rep("Control", nrow(ed_pat_c)),
                       rep("Pacientes", nrow(ed_pat_pat)))))
lims = c(min(toplot[, c("cronologica", "biologica")], na.rm = T), 
         max(toplot[, c("cronologica", "biologica")], na.rm = T))
ggplot(toplot, aes(x=cronologica, y=biologica, color = pacientes)) + 
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = c("Control" = "#31a7f6",
                                "Pacientes" = "#EB9050")) +
  # scale_color_manual(values = c("blue", "#FF7F50")) +
  xlab("Chronological Age") +
  ylab("Biological Age") +
  geom_abline(intercept = 0, linetype = 2, color = "darkgrey")+
  default_theme()+ 
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right") +
  coord_equal(xlim = lims, ylim = lims)


# Save -------------------------------------------------------------------------

save(bioage_c, bioage_pat, filt_bioage_c, filt_bioage_pat, 
     file = paste0(results_folder, "bioage.RData"))

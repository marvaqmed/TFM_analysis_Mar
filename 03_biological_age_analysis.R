# DESCRIPCIÓN: 
# Estimar la edad biológica del grupo control y de pacientes FIV
# Representación de la edad biológica frente a edad cronológica 
# Cálculo aceleraciones y diferencia de medias 
# Cálculo recta patrón en base a controles
# Estratificación de pacientes según aceleración epigenética relativa

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
source("Funciones_R.R")

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

ggsave(paste0(results_folder, "plots_controles.jpg"), dpi = 500)

#### Análisis de correlación ####
corrs = do.call("rbind", lapply(colnames(filt_bioage_c), function(x){
  res = cor.test(ed_pat_c$AGE2, filt_bioage_c[, x])
  data.frame(correlation = res$estimate, 
             p_value = res$p.value)
}))
rownames(corrs) = colnames(filt_bioage_c)
corrs
#          correlation      p_value
# Horvath    0.7646908 1.372159e-04
# BLUP       0.8985316 1.745236e-07
# EN         0.9536538 2.704041e-10
# AltumAge   0.7892321 5.890655e-05

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

ggsave(paste0(results_folder, "plots_controles_corr.jpg"), dpi = 500)


#### Cálculo MAE ###
mae_c = do.call("rbind", lapply(c("Horvath", "BLUP", "EN", "AltumAge"), function(x){
  diff = mean(abs(bioage_c[, x] - bioage_c$age))
  data.frame(mae = diff)
}))
rownames(mae_c) = c("Horvath", "BLUP", "EN", "AltumAge")
mae_c
#                mae
# Horvath   6.094583
# BLUP     14.990120
# EN        1.751881
# AltumAge 10.353524

# Seleccionamos el reloj Zhang EN 

# Análisis edad biológica pacientes --------------------------------------------

# Aunque hemos seleccionado el reloj Zhang EN probamos el resto de relojes
# Seleccionamos solo las columnas de los relojes
filt_bioage_pat = bioage_pat[, c("Horvath", "BLUP", "EN", "AltumAge")]

#### Edad biológica vs edad cronológica #### 
plots = lapply(colnames(filt_bioage_pat), function(x){
  toplot = data.frame(cronologica = ed_pat_pat$AGE2, 
                      biologica = filt_bioage_pat[, x])
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

ggsave(paste0(results_folder, "plots_pacientes.jpg"), dpi = 500)

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
  toplot = data.frame(cronologica = ed_pat_pat$AGE2, 
                      biologica = filt_bioage_pat[, x])
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

ggsave(paste0(results_folder, "plots_pacientes_corr.jpg"), dpi = 500)

#### Cálculo MAE ###
mae_pat = do.call("rbind", lapply(c("Horvath", "BLUP", "EN", "AltumAge"), function(x){
  diff = mean(abs(bioage_pat[, x] - bioage_pat$age))
  data.frame(mae = diff)
}))
rownames(mae_pat) = c("Horvath", "BLUP", "EN", "AltumAge")
mae_pat
#                mae
# Horvath   5.754229
# BLUP     15.778835
# EN        2.520933
# AltumAge  7.719855


# Pactientes FIV frente a controles --------------------------------------------

toplot = data.frame(
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
  xlab("Chronological Age") +
  ylab("Biological Age") +
  geom_abline(intercept = 0, linetype = 2, color = "darkgrey")+
  default_theme()+ 
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right") +
  coord_equal(xlim = lims, ylim = lims)

ggsave(paste0(results_folder, "cruce_rectas.jpg"), dpi = 500)

# Edad a la que se cruzan las rectas de controles y pacientes
# Recta controles
lm(formula = EN ~ age, data = bioage_c)
# y = 4.488 + 0.907x

# Recta pacientes 
lm(formula = EN ~ age, data = bioage_pat)
# y = 15.232 + 0.639x

# Punto de corte --> 40,22 años
# y = 4.488 + 0.907x = 15.232 + 0.6399x --> 0.907x - 0.6399x = 15.232 - 4.488 --> 
# 0.2671x = 10.744 --> x = 10.744 / 0.2671 = 40.22  



#### Diferencia media cronológica y biológica ####

# Cronológica controles
mean(bioage_c$age)
# 36.61587

# Biológica controles
mean(bioage_c$EN)
# 37.69823

# Cronológica pacientes
mean(bioage_pat$age)
# 40.38186

# Biológica pacientes
mean(bioage_pat$EN)
# 41.07115

# Diferencia medias controles
mean(bioage_c$EN) - mean(bioage_c$age)
# 1.082358

# Diferencia medias pacientes 
mean(bioage_pat$EN) - mean(bioage_pat$age)
# 0.6892878

# Determinar significancia estadística 
# Comprobamos si siguen distribución normal 
shapiro.test(bioage_c$age)
# p-value = 0.2261
shapiro.test(bioage_c$EN)
# p-value = 0.2116
shapiro.test(bioage_pat$age)
# p-value = 0.2486
shapiro.test(bioage_pat$EN)
# p-value = 0.6237

# Como todos siguen una distribución normal usamos t.test para comparar si 
# la edad biológica es significativamente mayor que la cronológica
t.test(bioage_c$age, bioage_c$EN, paired = T)
# p-value = 0.03138
t.test(bioage_pat$age, bioage_pat$EN, paired = T)
# p-value = 0.007114


#### Separando a los 40 años ####
# Cronológica controles -40.22
mean(bioage_c$age[bioage_c$age < 40.22])
# 33.27837

# Biológica controles -40.22
mean(bioage_c$EN[bioage_c$age < 40.22])
# 34.94852

# Cronológica controles +40.22
mean(bioage_c$age[bioage_c$age > 40.22])
# 43.84714

# Biológica controles +40.22
mean(bioage_c$EN[bioage_c$age > 40.22])
# 43.65594

# Diferencia medias controles -40.22
mean(bioage_c$EN[bioage_c$age < 40.22]) - mean(bioage_c$age[bioage_c$age < 40.22])
# 1.670153

# Diferencia medias controles +40.22
mean(bioage_c$EN[bioage_c$age > 40.22]) - mean(bioage_c$age[bioage_c$age > 40.22])
# -0.1911973

t.test(bioage_c$EN[bioage_c$age < 40.22] - bioage_c$age[bioage_c$age < 40.22], 
       bioage_c$EN[bioage_c$age > 40.22] - bioage_c$age[bioage_c$age > 40.22])
# p-value = 0.07744


## Pacientes 
# Cronológica pacientes -40.22
mean(bioage_pat$age[bioage_pat$age < 40.22])
# 36.51837

# Biológica pacientes -40.22
mean(bioage_pat$EN[bioage_pat$age < 40.22])
# 38.70515

# Cronológica pacientes +40.22
mean(bioage_pat$age[bioage_pat$age > 40.22])
# 43.67298

# Biológica pacientes +40.22
mean(bioage_pat$EN[bioage_pat$age > 40.22])
# 43.08662

# Diferencia medias pacientes -40.22
mean(bioage_pat$EN[bioage_pat$age < 40.22]) - mean(bioage_pat$age[bioage_pat$age < 40.22])
# 2.186788

# Diferencia medias pacientes +40
mean(bioage_pat$EN[bioage_pat$age > 40.22]) - mean(bioage_pat$age[bioage_pat$age > 40.22])
# -0.5863608

t.test(bioage_pat$EN[bioage_pat$age < 40.22] - bioage_pat$age[bioage_pat$age < 40.22], 
       bioage_pat$EN[bioage_pat$age > 40.22] - bioage_pat$age[bioage_pat$age > 40.22])
# p-value = 7.816e-09



# Estratificación pacientes -----------------------------------------------

#### Ajuste modelo con controles EN #### 
# lm controles
lm(formula = EN ~ age, data = bioage_c)
# interecpt = 4.488
# pendiente = 0.907
# y = x*0.907 + 4.488

# Para cada paciente sacar valor edad control 
# Y la aceleración respecto a la edad control
for (i in 1:nrow(bioage_pat)) {
  bioage_pat$control[i] = 0.907*(bioage_pat$age[i]) + 4.488
  bioage_pat$Acc2[i] = bioage_pat$EN[i] - bioage_pat$control[i]
}

# Estratificamos la población considerando un margen de error de +-1 año
bioage_pat$g_acc = ifelse(bioage_pat$EN > bioage_pat$control+1, "E", 
                        ifelse(bioage_pat$EN < bioage_pat$control-1, "R", "Normal"))
table(bioage_pat$g_acc)
#  E   Normal    R 
# 53       43   54

# Añadimos la estratificación a la metadata 
all(rownames(ed_pat_pat) == rownames(bioage_pat))
# TRUE 
ed_pat_pat$g_acc = bioage_pat$g_acc

# Representación gráfica
toplot = data.frame(cronologica = bioage_pat$age, 
                    biologica = bioage_pat$EN)
lims = c(min(toplot, na.rm = T), max(toplot, na.rm = T))
toplot$accel = factor(bioage_pat$g_acc, levels = c("E", "R", "Normal"))
ggplot(toplot, aes(x=cronologica, y=biologica, color = accel))+
  geom_point(size = 2) +
  scale_color_manual(values = c("#B75180", "#3FBCC3", "grey")) + 
  xlab("Edad cronológica")+
  ylab("Edad biológica")+
  ggtitle("Zhang EN")+
  geom_abline(intercept = 4.488, slope = 0.907, linetype = 2, color = "grey25") +
  geom_abline(intercept = 4.488 + 1, slope = 0.907, linetype = 3, color = "grey25") +
  geom_abline(intercept = 4.488 - 1, slope = 0.907, linetype = 3, color = "grey25") +
  default_theme() +
  theme(plot.title = element_text(hjust=0.5), 
        legend.position = "right", 
        axis.title = element_text(size = 20)) + 
  coord_equal(xlim = lims, ylim = lims)

ggsave("/data/local/mvaquerizo/metilacion/Plots/E_vs_R_total.jpg", dpi = 500,
       width = 10, height = 8)

# Save -------------------------------------------------------------------------

save(bioage_c, bioage_pat, ed_pat_pat, 
     file = paste0(results_folder, "bioage.RData"))

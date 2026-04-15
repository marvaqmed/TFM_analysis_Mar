# DESCRIPCIÓN: 
# Estudiar resultados reproductivos en función de la edad y la aceleración 
# en pacientes FIV

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
library(dplyr)

# Funciones 
source("../TFM/Funciones.R")

# Carpeta de resultados
results_folder = paste0(root, "results/04_reproductive_outcome/")
dir.create(results_folder)


# Carga de datos ---------------------------------------------------------------

load(paste0(root, "results/02_exploratory_analysis/data_filtered.RData"), verbose =T)
# Loading objects:
#   Betas
#   ed_exp
#   ed_pat
load(paste0(root, "results/03_biological_age_analysis/bioage.RData"), verbose = T)
# Loading objects:
#  bioage_c
#  bioage_pat
#  filt_bioage_c
#  filt_bioage_pat

#### Separar datos controles y pacientes ####
control = rownames(ed_pat[ed_pat$STUDY == 'REPROCYCLE', ])
length(control)
# 19

# Seleccionar solo la info de controles
betas_c = Betas[, control]
ed_pat_c = ed_pat[control, ]
ed_exp_c = ed_exp[control,]

# Seleccionar la info de las pacientes (!control)
betas_pat = Betas[, !colnames(Betas)%in%control]
ed_pat_pat = ed_pat[!rownames(ed_pat)%in%control, ]
ed_exp_pat = ed_exp[!rownames(ed_exp)%in%control, ]


# Filtrado de datos ------------------------------------------------------------
# Número de embriones transferidos
table(ed_pat_pat$TR_No_EMBRYOS, useNA = "ifany")
#   1  2 NA 
# 143  1  6
# Hay uno que no es SET -> PAT_8
# Lo eliminamos 
ed2_pat_pat = ed_pat_pat[rownames(ed_pat_pat)!="PAT_8", ]
dim(ed2_pat_pat)
# 149 35

# Día del embrión
table(ed2_pat_pat$TR_EMBRYO_DAY)
#  D5  D6
# 120  29
# Nos quedamos solo con D5 
ed2_pat_pat = ed2_pat_pat[ed2_pat_pat$TR_EMBRYO_DAY =="D5", ]
dim(ed2_pat_pat)
# 120  35

# Calidad del embrión
table(ed2_pat_pat$TR_EMBRYO_QUALITY)
#  A  B  C 
# 23 89  8
# Nos quedamos con embriones A y B
ed2_pat_pat = ed2_pat_pat[ed2_pat_pat$TR_EMBRYO_QUALITY != "C", ]
dim(ed2_pat_pat)
# 112  35

# Origen ovocito
table(ed2_pat_pat$TR_OOCYTE_ORIGIN, ed2_pat_pat$TR_PGT)
#       Euploid  NO_PGT
# Donated     7      59
# Own        41       5
table(ed2_pat_pat$TR_PGT, ed2_pat_pat$AGE)
#         32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
# Euploid  1  3  4  3  2  6  4  5  8  5  2  1  4  0  1  0  0  0
# NO_PGT   3  1  1  3  1  1  1  7  3  8  2  4  7  6  6  5  2  2
# Nos quedamos con ovodonación, propio + PGT o propio <35 años
ed2_pat_pat = ed2_pat_pat[(ed2_pat_pat$TR_OOCYTE_ORIGIN == "Donated") | 
                            (ed2_pat_pat$TR_OOCYTE_ORIGIN == "Own" & ed2_pat_pat$TR_PGT =="Euploid") |
                            (ed2_pat_pat$TR_OOCYTE_ORIGIN == "Own" & ed2_pat_pat$AGE2 < 35), ]
dim(ed2_pat_pat)
# 110  35

# Juntamos miscarriages en un único level
ed2_pat_pat$TR_OUTCOME = as.character(ed2_pat_pat$TR_OUTCOME)
ed2_pat_pat$TR_OUTCOME[ed2_pat_pat$TR_OUTCOME %in% c("BM", "CM")] = "M"
ed2_pat_pat$TR_OUTCOME = factor(ed2_pat_pat$TR_OUTCOME, levels = c("LB", "M", "NB"))
table(ed2_pat_pat$TR_EMBRYO_QUALITY, ed2_pat_pat$TR_OUTCOME)
#   LB  M NB
# A 11  5  7
# B 43 14 30

# Seleccionamos las mismas filas en el df de edad
filt2_bioage_pat = filt_bioage_pat[rownames(ed2_pat_pat), ]
dim(filt2_bioage_pat)
# 110 4
bioage2_pat = bioage_pat[rownames(ed2_pat_pat), ]
dim(bioage2_pat)
# 110 12


# Exploración  -----------------------------------------------------------------

#### Correlación controles y pacientes #### 
# Representación gráfica
toplot <- data.frame(
  cronologica = c(ed_pat_c$AGE2, ed2_pat_pat$AGE2),
  biologica = c(filt_bioage_c$EN, filt2_bioage_pat$EN),
  pacientes = factor(c(rep("Control", nrow(ed_pat_c)),
                       rep("Pacientes", nrow(ed2_pat_pat)))))
lims = c(min(toplot[, c("cronologica", "biologica")], na.rm = T), max(toplot[, c("cronologica", "biologica")], na.rm = T))
ggplot(toplot, aes(x=cronologica, y=biologica, color = pacientes)) + 
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = c("blue", "#FF7F50")) +
  xlab("Chronological Age") +
  ylab("Biological Age") +
  geom_abline(intercept = 0, linetype = 2, color = "darkgrey")+
  default_theme()+ 
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right") +
  coord_equal(xlim = lims, ylim = lims)

# Correlación 
cor.test(ed_pat_c$AGE2, filt_bioage_c$EN)
# 0.9537
cor.test(ed2_pat_pat$AGE2, filt2_bioage_pat$EN)
# 0.676


#### Ajuste modelo con controles EN #### 
datos_EN = cbind(filt2_bioage_pat[, c("EN")], ed2_pat_pat)
colnames(datos_EN)[1] = "EN"
# lm controles
lm(formula = EN ~ age, data = bioage_c)
# interecpt = 4.488
# pendiente = 0.907
# y = x*0.907 + 4.488

# Para cada paciente sacar valor edad control 
# Y la aceleración respecto a la edad control
for (i in 1:nrow(datos_EN)) {
  datos_EN$control[i] = 0.907*(datos_EN$AGE2[i]) + 4.488
  datos_EN$Acc2[i] = datos_EN$EN[i] - datos_EN$control[i]
}

# Estratificamos la población considerando un margen de error de +-1 año
datos_EN$g_acc = ifelse(datos_EN$EN > datos_EN$control+1, "E", 
                        ifelse(datos_EN$EN < datos_EN$control-1, "R", "Normal"))
table(datos_EN$g_acc)
#  E++   Normal    R++ 
#   39       29     42 

# Representación gráfica
toplot = data.frame(cronologica = datos_EN$AGE2, 
                    biologica = datos_EN$EN)
lims = c(min(toplot, na.rm = T), max(toplot, na.rm = T))
toplot$accel = factor(datos_EN$g_acc, levels = c("E", "R", "Normal"))
ggplot(toplot, aes(x=cronologica, y=biologica, color = accel))+
  geom_point(size = 2) +
  scale_color_manual(values = c("#B75180", "#3FBCC3", "grey")) + 
  xlab("Edad cronológica")+
  ylab("Edad biológica (Zhang EN)")+
  ggtitle("")+
  geom_abline(intercept = 4.488, slope = 0.907, linetype = 2, color = "grey25") +
  geom_abline(intercept = 4.488 + 1, slope = 0.907, linetype = 3, color = "grey25") +
  geom_abline(intercept = 4.488 - 1, slope = 0.907, linetype = 3, color = "grey25") +
  default_theme() +
  theme(plot.title = element_text(hjust=0.5), 
        legend.position = "right", 
        axis.title = element_text(size = 20)) + 
  coord_equal(xlim = lims, ylim = lims)


#### Tabla demográfica #### 
demo_dat = ed2_pat_pat[, c("AGE2", "BMI", "BUP", "TR_EMBRYO_STATE", "TR_EMBRYO_QUALITY", 
                           "TR_OOCYTE_ORIGIN", "TR_OOCYTE_STATE", "TR_CYCLE_TYPE", 
                           "TR_PGT")]
demo_dat$TR_EMBRYO_STATE = ifelse(demo_dat$TR_EMBRYO_STATE == "Fresh", "Fresh", "Frozen")
demo_dat$TR_OOCYTE_STATE = factor(ifelse(demo_dat$TR_OOCYTE_STATE == "Fresh", "Fresh", "Frozen"), 
                                  levels = c("Fresh", "Frozen"))
demo_dat$TR_CYCLE_TYPE = factor(ifelse(demo_dat$TR_CYCLE_TYPE == "Natural", "Natural", "HRT"), 
                                levels = c("Natural", "HRT"))
demo_dat$g_acc = datos_EN$g_acc
tabla = demographic_table(demo_dat, demo_dat$g_acc)
#write.csv(tabla, file = paste0(results_folder, "tabla_demografica_E_vs_R.csv"))




# Asociación con resultados reproductivos --------------------------------------

# Seleccionamos solo la información de envejecidas y rejuvenecidas
# NO tenemos en cuenta las pacientes cuya aceleración está en el rango
sel_datos = datos_EN[datos_EN$g_acc != "Normal", ]

#### Tasas reproductivas ####
tasas = rbind(LBR1 = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] == "LB")/sum(sel_datos$g_acc == "E"),
                                Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] == "LB")/sum(sel_datos$g_acc == "R")), 
              LBR2 = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] == "LB")/sum(sel_datos$outcome[sel_datos$g_acc == "E"] %in% c("LB", "M")),
                                Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] == "LB")/sum(sel_datos$outcome[sel_datos$g_acc == "R"] %in% c("LB", "M"))),
              MR1 = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] == "M")/sum(sel_datos$g_acc == "E"), 
                               Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] == "M")/sum(sel_datos$g_acc == "R")), 
              MR2 = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] == "M")/sum(sel_datos$outcome[sel_datos$g_acc == "E"] %in% c("LB", "M")), 
                              Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] == "M")/sum(sel_datos$outcome[sel_datos$g_acc == "R"] %in% c("LB", "M"))), 
              NB = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] == "NB")/sum(sel_datos$g_acc == "E"), 
                              Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] == "NB")/sum(sel_datos$g_acc == "R")), 
              IR = data.frame(Envejecida = sum(sel_datos$outcome[sel_datos$g_acc == "E"] %in% c("LB", "M"))/sum(sel_datos$g_acc == "E"), 
                              Rejuvenecida = sum(sel_datos$outcome[sel_datos$g_acc == "R"] %in% c("LB", "M"))/sum(sel_datos$g_acc == "R"))
)
tasas

#      Envejecida Rejuvenecida
# LBR1  0.4615385    0.5238095
# LBR2  0.7200000    0.7586207
# MR1   0.1794872    0.1666667
# MR2   0.2800000    0.2413793
# NB    0.3589744    0.3095238
# IR    0.6410256    0.6904762

# P-valores
LBR1 = ifelse(sel_datos$outcome == "LB", "LB", "noLB")
MR1 = ifelse(sel_datos$outcome == "M", "M", "noM")
NB = ifelse(sel_datos$outcome == "NB", "NB", "noNB")
IR = ifelse(sel_datos$outcome %in% c("LB", "M"), "I", "noI")

grupos = list(LBR = LBR1, MR = MR1, NB = NB, IR = IR)

pvals = do.call("rbind", lapply(grupos, function(x){
  res = fisher.test(table(x, sel_datos$g_acc))
  data.frame(p_value = res$p.value)
}))
pvals
#       p_value
# LBR 0.6586002
# MR  1.0000000
# NB  0.6465124
# IR  0.6465124


# Tasas apiladas 
sel_datos$g_acc = factor(sel_datos$g_acc, levels = c('R', 'E'))
aa = sweep(as.matrix(table(sel_datos$TR_OUTCOME, sel_datos$g_acc)), 2, 
           as.numeric(table(sel_datos$g_acc)), FUN ='/')

toplot2 = melt(aa)
colnames(toplot2) = c("Outcome", "Acceleration", "value")
ggplot(toplot2, aes(x = Acceleration, y = value, fill = Outcome))+
  geom_bar(stat= "identity") +
  scale_fill_manual(values = c("#42BA97", "#E6B729", "#EA6312")) +
  default_theme()


pvals = do.call("rbind", lapply(grupos, function(x){
  res = fisher.test(table(x, datos_EN$group_accel))
  data.frame(p_value = res$p.value)
}))
pvals


# Save --------------------------------------------------------------------

save(datos_EN, sel_datos, file = paste0(results_folder, "datos_E_R.RData"))




default_theme = function(x.axis.angle = 0, legend.pos = "bottom"){
  
  hjust <- ifelse(x.axis.angle > 0, 1, 0)
  
  custom_theme  <- theme_light() +
    theme(legend.position = legend.pos,
          axis.title = element_text(size=18),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22, hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=14))
  
  if (abs(x.axis.angle) == 90){
    custom_theme <- custom_theme + theme(axis.text.x = element_text(angle=x.axis.angle, vjust = 0.5, hjust = hjust))
    } 
  else if (abs(x.axis.angle) == 45){
    custom_theme <- custom_theme + theme(axis.text.x = element_text(angle=x.axis.angle, vjust = 1, hjust = hjust))
    }
  
  return(custom_theme)
}

demographic_table = function(demog_dat, classification, levels = NULL) {
  
  # Fix data
  demog_dat[demog_dat == ''] = NA
  demog_dat[demog_dat == "NA"] = NA
  
  if (is.null(levels)) {classification = factor(classification)}
  else {classification = factor(classification, levels = levels)}
  
  # Unimos el data frame que lleva los datos demográficos que queremos junto con la
  # classification en el mismo data frame y quitamos los NAs
  demog_dat = cbind(demog_dat, classification = classification)
  demog_dat = demog_dat[!is.na(demog_dat$classification), ]
  
  # Vamos construyendo la tabla indicando lo que tiene que ir en cada columna
  # Añadimos la n y el porcentaje con respecto al total que supone cada grupo
  N_g = sum(!is.na(classification))
  
  res_g = cbind(
    N_g,
    t(paste0(
      table(classification),
      " (",
      round(t(table(classification)) * 100 / sum(!is.na(classification)), 1),
      "%)"
    )),
    "N/A"
  )
  colnames(res_g) = c("N", levels(classification), "p_value")
  
  # Generamos las filas con cada variable y hacemos el análisis demografico de cada una
  res = do.call(
    "rbind",
    lapply(colnames(demog_dat)[-ncol(demog_dat)], function(col) {
      cat(col, "\n")
      
      demog_dat_aux = demog_dat[!is.na(demog_dat[, col]), ]
      
      N = paste0(
        nrow(demog_dat_aux),
        " (",
        round(nrow(demog_dat_aux) * 100 / N_g, 1),
        "%)"
      ) # Columna N
      
      # Ahora definimos lo que va a ir en el resto de columnas dependiendo de si la
      # variable  es numérica (median +/- se e intervalo de confianza) o discreta
      # (la N de cada categoría). Dentro de la numérica, diferenciamos el tipo de test
      # utilizado dependiendo de si es una comparación simple o múltiple
      
      # Variables numéricas
      
      if (is.numeric(demog_dat_aux[, col])) {
        # Si la variable es numérica...
        # Definimos lo que queremos que vaya en las columnas de las estadísticas y el
        # número de decimales
        
        sume = summarySE(
          data = demog_dat_aux,
          groupvars = "classification",
          measurevar = col,
          na.rm = T
        )
        stats = apply(sume[, -1], 1, function(x) {
          paste0(
            round(x[col], 2),
            "+/-",
            round(x["sd"], 2),
            " [",
            round(x[col] - x["ci"], 2),
            ",",
            round(x[col] + x["ci"], 2),
            "]"
          )
        })
        names(stats) = levels(classification)
        
        # Si la variable tiene solo dos niveles, aplicamos un wilcox test e indicamos que
        # en la columna va el p-valor
        print(shapiro.test(demog_dat_aux[, col]))
        if (shapiro.test(demog_dat_aux[, col])$p.value > 0.05) {
          print(paste0("Normal:", col))
          
          if (length(unique(classification)) == 2) {
            tt = t.test(
              demog_dat_aux[
                demog_dat_aux$classification ==
                  levels(demog_dat_aux$classification)[1],
                col
              ],
              demog_dat_aux[
                demog_dat_aux$classification ==
                  levels(demog_dat_aux$classification)[2],
                col
              ]
            )
            p_value = tt$p.value
          } else {
            # Si tiene más de dos niveles, aplicamos un kruskal
            anv = aov(demog_dat[, col] ~ demog_dat$classification)
            p_value = summary(anv)[[1]][["Pr(>F)"]][[1]]
          }
        } else {
          if (length(unique(classification)) == 2) {
            tt = wilcox.test(
              demog_dat_aux[
                demog_dat_aux$classification ==
                  levels(demog_dat_aux$classification)[1],
                col
              ],
              demog_dat_aux[
                demog_dat_aux$classification ==
                  levels(demog_dat_aux$classification)[2],
                col
              ]
            )
            p_value = tt$p.value
          } else {
            # Si tiene más de dos niveles, aplicamos un kruskal
            kw = kruskal.test(demog_dat[, col], demog_dat$classification)
            p_value = kw$p.value
          }
        }
      } else {
       
        # Variables discretas
        # Si la variable es discreta, le indicamos que ponga las N de cada grupo separados
        # por ;
        
        tt = table(demog_dat_aux[, col], demog_dat_aux$classification)
        aux2 = do.call(
          "rbind",
          lapply(rownames(tt), function(x) {
            aux = paste0(
              x,
              " = ",
              tt[x, ],
              " (",
              round(tt[x, ] * 100 / table(demog_dat_aux$classification), 1),
              "%)"
            )
          })
        )
        stats = apply(aux2, 2, paste, collapse = "; ")
        names(stats) = levels(classification)
        if (nrow(tt) == 1) {
          # Para la columna del p-valor, le decimos que en el caso de que
          # todas las pacientes tengan el mismo valor para esa variable (p.ej. ciclo natural),
          # tiene que poner N/A y no hacer el test
          p_value = "N/A"
        } else {
          # En caso de que tenga más de un grupo, le decimos que haga un fisher test
          ft = fisher.test(tt, simulate.p.value = T) # Si hay alguna columna grande que cueste mucho hacer poner, simulate.p.value=TRUE
          p_value = ft$p.value
        }
      }
      
      result = c(N = N, stats, p_value = p_value) # Indicamos lo que tiene que ir en cada columna
      return(result)
    })
  )
  
  rownames(res) = colnames(demog_dat)[-ncol(demog_dat)]
  
  res = rbind(res_g, res)
  return(res)
}



plot_PCAscores = function(dat, ed, components = c(1,2), condition1, colors = NULL, condition2=NULL, shapes=NULL,
                          title = "", legendpos = "bottom", labels = NULL, ellipses = F, ellipse_level=0.95,
                          ellipse_alpha=0.2){

  ## Ejecutamos el PCA
  pca = prcomp(t(dat))
  Var = round(summary(pca)$importance[2,components]*100,1)
  
  ## Generamos data frame con las puntuaciones para representarlo
  toplot = data.frame(pca$x[,components], stringsAsFactors = F)
  # Añadimos color en base a la condition1
  toplot$color = ed[colnames(dat), condition1]
  # Paleta por defecto
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7", "#9370DB", "#B3EE3A", "#FFC1C1",
                 "#999999", "#000000")
  if(is.null(colors)){
    colors = cbPalette[1:length(unique(ed[, condition1]))]
  }
  # Añadimos forma de los puntos en base a condition2
  if (!is.null(condition2)){
    toplot$shape = ed[colnames(dat), condition2]
    # Formas por defecto
    if(is.null(shapes)){
      shapes = 15:(length(unique(ed[, condition2]))+15)
    }
  }
  # Añadimos nombre de las muestras
  if (!is.null(labels)){
    toplot$labels = ed[colnames(dat), labels]
  }
  
  ## Límites de los ejes
  lim = max(abs(c(min(toplot), max(toplot))))
  # Si ellipses = T los límites los hacemos más grandes
  if (ellipses){
    lim = lim + 0.8*lim
  }
  axis_limits = c(-lim, lim)
  
  ## Generamos el plot 
  p = ggplot(toplot, aes_string(x=colnames(toplot)[1], y=colnames(toplot)[2]))
  # Colores y formas
  if (!is.null(condition2)){
    p = p + geom_point(aes(color=color, shape=shape), size=3) +
      scale_color_manual(name = condition1, values = colors)+
      scale_shape_manual(name = condition2, values = shapes)
  }else{
    p = p + geom_point(aes(color=color), size=3) +
      scale_color_manual(name = condition1, values = colors)
  }
  # Nombre de las muestras
  if(!is.null(labels)){
    p = p + geom_text_repel(aes(label=labels, color=toplot$color), size=4,
                            point.padding = unit(0.2, 'lines'),segment.size = 0 )
  }
  
  if (ellipses){
    p = p + stat_ellipse(aes(fill=color),geom="polygon",level=ellipse_level, alpha=ellipse_alpha,
                         show.legend = FALSE) +
      scale_fill_manual(values=colors)
  }
  # Información de los ejes y título
  p = p + xlab(paste0("PC", components[1], ": ", Var[1], "%"))+
    ylab(paste0("PC", components[2], ": ", Var[2], "%"))+
    ggtitle(title)
  # Límites de los ejes
  p = p + xlim(axis_limits) + ylim(axis_limits)
  # Parámetros gráficos
  p = p + theme_light() +
    theme(legend.position = legendpos,
          axis.title = element_text(size=18),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22, hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=13))
  return(p)
}


plot_Boxplot = function(dat, ed, condition = NULL, colors = NULL, title = "", legendpos = "bottom"){
  
  # Paleta de colores por defecto
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7", "#9370DB", "#B3EE3A", "#FFC1C1",
                 "#999999", "#000000")
  if(is.null(colors)){
    colors = cbPalette[1:length(unique(ed[, condition]))]
  }
  
  # Si no hay condición, se pintan todas las barras grises
  if(is.null(condition)){
    colors = rep("#999999", ncol(dat))
  }
  
  # Plot boxplots
  toplot = melt(as.matrix(dat))
  colnames(toplot) = c("gene", "variable", "value")
  toplot$group = ed[toplot$variable, condition]
  p = ggplot(toplot, aes(x=variable, y=value, fill=group))+
    geom_boxplot()+
    scale_fill_manual(name = condition, values = colors)+
    xlab("")+ ylab("")+ ggtitle(title)+
    theme_light()+
    theme(legend.position = legendpos,
          axis.title = element_text(size=18),
          axis.text = element_text(size=15),
          axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
          plot.title = element_text(size=22, hjust=0.5),
          legend.title = element_blank(),
          legend.text = element_text(size=13))
  return(p)
}


diffExprAnalysis = function(dat, ed, condition, paired = FALSE, pair = NULL){
  
  # Creamos experimental design
  group = ed[colnames(dat), condition]
  if (paired){
    pair = ed[colnames(dat), pair]
    design <- model.matrix(~ 0+group+pair)
    rownames(design) <- colnames(dat)
    fit <- lmFit (dat, design= design)
    res.limma <- eBayes (fit)
    res = list(topTable(res.limma, coef=1, number = nrow((dat))))
    names(res) = paste(gsub("group", "", colnames(design)[1:2]), collapse="-")
  }else{
    design <- model.matrix(~ 0+group)
    rownames(design) <- colnames(dat)
    
    # Generamos los contrastes
    combs = combn(x=colnames(design), m=2, simplify = TRUE)
    cons2 = apply(combs,2,function(y){
      paste(y,collapse="-")
    })
    cons = gsub("group","",cons2)
    contrasts <- makeContrasts(contrasts=cons2, levels = design)
    
    # Ajuste modelo lineal
    fit <- lmFit (dat, design= design)
    fit.cont <- contrasts.fit (fit, contrasts)
    res.limma <- eBayes (fit.cont)
    
    # Extraemos los resultados
    res = list()
    for(i in 1:length(cons2)){
      res[[cons2[i]]] = topTable(res.limma, coef = i, number = nrow(dat))
    }
    names(res) = gsub("group","",names(res))
  }
  
  # Incluimos columna FDRs 
  res2 = lapply(names(res), function(x){
    aux = unlist(strsplit(x, split = "-", fixed = T))
    condition1 = aux[1]
    condition2 = aux[2]
    res[[x]]$FC = unlist(lapply(rownames(res[[x]]),function(y){
      A = mean(as.numeric(dat[y, group==condition1]))
      B = mean(as.numeric(dat[y, group==condition2]))
      FC = 2^(A-B)
      if (FC<1){
        FC = -1/FC
      }
      FC
    }))
    res[[x]]
  })
  names(res2) = names(res)
  
  return(res2)
}
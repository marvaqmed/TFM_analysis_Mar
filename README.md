# Edad molecular endometrial en infertilidad: relojes epgenéticos y consecuencias funcionales 
## Descripción del proyecto 
Este repositorio contiene el flujo de trabajo bioinformático y estadístico desarrollado para el Trabajo Fin de Máster (TFM) del Máster en Bioinformática de la Universitat de València. 

El objetivo principal del estudio es determinar cuál es el reloj epigenético más adecuado para evaluar el envejecimiento endometrial en el contexto de la reproducción asistida. 

Todos los análisis se han realizado combinando los entornos de programación R y Python. El flujo de análisis se inicia a partir de datos de metilación del ADN obtenidos de biopsias endometriales procedentes tanto por pacientes sometidas a tratamientos de fecundación in vitro como de voluntarias sin sospecha de infertilidad incluidas como grupo control. 
Tras el procesamiento y control de calidad de los datos de metilación, se estimó la edad biológica endometrial mediante distintos relojes epigenéticos. En concreto, los relojes de Horvath, Zhang EN y Zhang BLUP fueron implementados en R, mientras que el reloj AltumAge fue aplicado en Python. 
Posteriormente, la cohorte de pacientes fue estratificada en base a la aceleración epigenética relativa y los grupos obtenidos fueron caracterizados clínica y funcionalmente. Para esto último, se integraron datos transcriptómicos obtenidos mediante RNA-seq, los cuales fueron procesados y normalizados previamente.
A continuación, se llevó a cabo un análisis de expresión diferencial utilizando el paquete limma, seguido de un análisis de enriquecimiento de conjunto de genes utilizando la anotación funcional de las bases de datos GeneOntology (GO) y Kyoto Encyclopedia of Genes and Genomes (KEGG). 
Finalmente, se aplicaron modelos de modelos de aprendizaje automático en Python para la identificación de biomarcadores asociados al envejecimiento epigenético endometrial. 

--- 
## Organización del repositorio 
```text
.
├── Funciones_python/                  
│   └── Funciones auxiliares utilizadas en Python.
│
├── pyaging_data/                      
│   └── Archivos necesarios para los análisis con pyaging.
│
├── 01_fromRAW_to_beta.R               
│   └── Preprocesamiento de los datos de metilación para la obtención de la matriz de betas.
│
├── 02_exploratory_analysis.R          
│   └── Análisis exploratorio y control de calidad de los datos de metilación.  
│
├── 03.1_AltumAge.ipynb                
│   └── Estimación de la edad biológica mediante AltumAge.
│
├── 03_biological_age_analysis.R       
│   └── Estimación de la edad epigenética mediante methylclock y estratificación de la población en base a la aceleración epigenética relativa.
│
├── 04_reproductive_outcomes.R         
│   └── Caracterización clínica de los grupos de envejecimiento epigenético y su asociación con los resultados reproductivos.
│
├── 05_select_transcriptome_samples.R  
│   └── Selección de muestras con datos de transcriptómica.
│
├── 06_exploratory_and_DEA_limma.R     
│   └── Análisis exploratorio de los datos transcriptómicos y análisis de expresión diferencial mediante limma.
│
├── 07_functional_analysis.R           
│   └── Análisis funcional mediante GSEA y anotación con GO y KEGG.
│
├── 08_gene_signature.ipynb            
│   └── Identificación de biomarcadores asociados a envejecimiento epigenético endometrial mediante la aplicación de modelos de aprendizaje automático.
│
├── Funciones_R.R                      
│   └── Funciones auxiliares utilizadas en R.
│
└── README.md
```

## Requisitos 

### Entorno de trabajo 
- R v.4.4.1
- Python v.3.12.3

### Principales paquetes de R 
- minfi 	
- ewastools 
- methylclock
- pvca
- HotellingEllipse
- FactoMineR
- VGAM
- limma
- NOISeq
- clusterProfiler
- org.Hs.eg.db
- ggplot2
- dplyr
- readxl
- ggvenn

### Principales paquetes de Python 
- numpy
- pandas
- pyaging 
- scikit-learn
- scipy
- igraph
- matplotlib
- seaborn

## Disponibilidad de los datos 
Los datos utilizados en este proyecto contienen información clínica y molecular sensible y, por lo tanto, no se encuentran disponibles públicamente debido a restricciones éticas y de confidencialidad. 
En este repositorio solo se proporciona el código necesario para reproducir los análisis bioinformáticos y estadísticos desarrollados en el estudio. 

## Uso 
Los scripts han sido diseñados para ejecutarse secuencialmente siguiendo el orden numérico del repositorio. 

## Autora 
Mar Vaquerizo Medina

Máster en Bioinformática 

Universitat de València 

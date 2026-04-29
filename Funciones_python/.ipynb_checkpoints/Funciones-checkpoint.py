# Load libraries

import pandas as pd
import statistics as st
import scipy.stats as sc
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix

# Funciones
def pca_plot(data, y:list):
    """
    Esta función dibuja un PCA indicando el % de variabilidad de cada componente
    Es específica para envejecidas y rejuvenecidas
    """
    x = StandardScaler().fit_transform(data)
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, 
                               columns=['PC1', 'PC2'])
    principalDf["label"] = y

    var_exp = pca.explained_variance_ratio_ * 100
    
    plt.figure(figsize=(10,10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    plt.xlabel(f'PC1 ({var_exp[0]:.2f}%)', fontsize=20)
    plt.ylabel(f'PC2 ({var_exp[1]:.2f}%)', fontsize=20)
    plt.title("Principal Component Analysis", fontsize=20)
    
    targets = ["E", "R"]
    colors = ["#3FBCC3", "#B75180"]
    
    for target, color in zip(targets, colors):
        indicesToKeep = principalDf['label'] == target
        plt.scatter(principalDf.loc[indicesToKeep, 'PC1'],
                    principalDf.loc[indicesToKeep, 'PC2'],
                    c = color, 
                    s = 50)
    
    plt.legend(targets, prop={'size': 15})
    plt.grid()
    plt.show()

def corr_loop(mytupleindex):
    """
    Esta función calcula los coeficientes de correlación entre pares de genes
    """
    subdat, i, j = mytupleindex
    res = sc.pearsonr(subdat.iloc[i,:], subdat.iloc[j,:])
    res2 = pd.DataFrame.from_dict({'Gene1': [subdat.index[i]], 'Gene2':[subdat.index[j]],  
                                'R': [res.statistic], 'p_value': [res.pvalue]})
    return(res2)


def accuracy_models(cross_val, data:pd.DataFrame, y:list): 
    """
    Esta función lanza y calcula el accuracy medio de los modelos de clasificación: SVM, Random Forest y KNN
    """
    accuracy_svm = list()
    accuracy_rf = list()
    accuracy_knn = list()
    
    for time, value in cross_val.items():
        
        for fold, val2 in value.items():

            # SVM 
            model_svm = SVC(random_state=123, probability=True)
            model_svm.fit(data.iloc[val2['Train'],:], np.array(y)[val2['Train']])
            y_pred = np.argmax(model_svm.predict_proba(data.iloc[val2['Test'],:]), axis = 1)
            y_true = [1 if label == 'R' else 0 for label in np.array(y)[val2['Test']]]
            tn, fp, fn, tp = confusion_matrix(y_pred,y_true).ravel().tolist()
            accuracy_svm.append((tn+tp)/(tn+tp+fp+fn))
                            
            # Random Forest
            model_rf = RandomForestClassifier(random_state=123)
            model_rf.fit(data.iloc[val2['Train'],:], np.array(y)[val2['Train']])
            y_pred = np.argmax(model_rf.predict_proba(data.iloc[val2['Test'],:]), axis = 1)
            y_true = [1 if label == 'R' else 0 for label in np.array(y)[val2['Test']]]
            tn, fp, fn, tp = confusion_matrix(y_pred,y_true).ravel().tolist()
            accuracy_rf.append((tn+tp)/(tn+tp+fp+fn))

            # KNN 
            model_knn = KNeighborsClassifier()
            model_knn.fit(data.iloc[val2['Train'],:], np.array(y)[val2['Train']])
            y_pred = np.argmax(model_knn.predict_proba(data.iloc[val2['Test'],:]), axis = 1)
            y_true = [1 if label == 'R' else 0 for label in np.array(y)[val2['Test']]]
            tn, fp, fn, tp = confusion_matrix(y_pred,y_true).ravel().tolist()
            accuracy_knn.append((tn+tp)/(tn+tp+fp+fn))
            

    return st.mean(accuracy_svm), st.mean(accuracy_rf), st.mean(accuracy_knn)

def prob_to_clase(probs:list, umbral:0.5): 
    """
    Saca la clase más probable (E - R) a partir de la lista de probabilidades
    """
    return ['E' if p > umbral else 'R' for p in probs]

def performance_modelo(y_true:list, y_pred:list, performance, modelo:str): 
    tn, fp, fn, tp = confusion_matrix(y_pred, y_true).ravel().tolist()
    performance.loc['Acc', modelo] = round((tn+tp)/(tn+tp+fp+fn), 3)
    performance.loc['Sens', modelo] = round((tp)/(tp+fn), 3)
    performance.loc['Spec', modelo] = round(tn/(tn+fp), 3)
    performance.loc['F1', modelo] = round((2*tp)/(2*tp+fp+fn), 3)
    


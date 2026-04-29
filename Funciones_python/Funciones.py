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
    colors = ["#B75180", "#3FBCC3"]
    
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
    Esta función lanza y calcula el accuracy medio de los modelos de clasificación: SVM, Random Forest y kNN
    """
    accuracy_svm = list()
    accuracy_rf = list()
    accuracy_knn = list()
    accuracy_svm_rf = list()
    accuracy_svm_knn = list()
    accuracy_rf_knn = list()
    
    for time, value in cross_val.items():
        
        for fold, val2 in value.items():

            # Separar datos
            data_train = data.iloc[val2['Train'],:]
            y_train = np.array(y)[val2['Train']]
            data_test = data.iloc[val2['Test'],:]
            y_test = [1 if label == 'R' else 0 for label in np.array(y)[val2['Test']]]

            # SVM 
            model_svm = SVC(random_state=123, probability=True)
            model_svm.fit(data_train, y_train)
            prob_svm = model_svm.predict_proba(data_test)
            y_pred = np.argmax(prob_svm, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_svm.append((tn+tp)/(tn+tp+fp+fn))
                            
            # Random Forest
            model_rf = RandomForestClassifier(random_state=123)
            model_rf.fit(data_train, y_train)
            prob_rf = model_rf.predict_proba(data.iloc[val2['Test'],:])
            y_pred = np.argmax(prob_rf, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_rf.append((tn+tp)/(tn+tp+fp+fn))

            # kNN 
            model_knn = KNeighborsClassifier()
            model_knn.fit(data_train, y_train)
            prob_knn = model_knn.predict_proba(data.iloc[val2['Test'],:])
            y_pred = np.argmax(prob_knn, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_knn.append((tn+tp)/(tn+tp+fp+fn))

            # SVM_Rf
            prob_svm_rf = np.mean([prob_svm, prob_rf], axis = 0)
            y_pred = np.argmax(prob_svm_rf, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_svm_rf.append((tn+tp)/(tn+tp+fp+fn))

            # SVM_kNN 
            prob_svm_knn = np.mean([prob_svm, prob_knn], axis = 0)
            y_pred = np.argmax(prob_svm_knn, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_svm_knn.append((tn+tp)/(tn+tp+fp+fn))   

            # Rf_kNN
            prob_rf_knn = np.mean([prob_rf, prob_knn], axis = 0)
            y_pred = np.argmax(prob_rf_knn, axis = 1)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel().tolist()
            accuracy_rf_knn.append((tn+tp)/(tn+tp+fp+fn))         

            

    return st.mean(accuracy_svm), st.mean(accuracy_rf), st.mean(accuracy_knn)


def prob_to_clase(probs:list): 
    """
    Saca la clase más probable (E - R) a partir de la lista de probabilidades
    """
    return ['E' if p >= 0.5 else 'R' for p in probs]


def performance_modelo(y_true:list, y_pred:list, performance, modelo:str): 
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=['R', 'E']).ravel().tolist()
    performance.loc['Acc', modelo] = round((tn+tp)/(tn+tp+fp+fn), 3)
    performance.loc['Sens', modelo] = round((tp)/(tp+fn), 3)
    performance.loc['Spec', modelo] = round(tn/(tn+fp), 3)
    performance.loc['F1', modelo] = round((2*tp)/(2*tp+fp+fn), 3)
    

def modelo_votos(probabilities: pd.DataFrame):
    votos = (probabilities[['SVM','Rf','kNN']] >= 0.5).sum(axis=1)
    
    return ['E' if voto >= 2 else 'R' for voto in votos]


def rendimiento_train(cross_val, data, y): 
    """
    Esta función lanza y calcula el rendimiento medio de los modelos de clasificación: SVM, Random Forest y kNN
    y de los combinados con los datos del entrenamiento.
    En cada time rellena el df de probabilidades con los 3 modelos y al final saca los modelos combinados y el 
    rendimiento de todos. 
    """

    accuracy = {'SVM': list(), 'Rf': list(), 'kNN': list(), 'SVM_Rf': list(), 'SVM_kNN': list(), 'Rf_kNN': list(), 'Prob_media': list(), 'Votos': list()}
    sensibilidad = {'SVM': list(), 'Rf': list(), 'kNN': list(), 'SVM_Rf': list(), 'SVM_kNN': list(), 'Rf_kNN': list(), 'Prob_media': list(), 'Votos': list()}
    especificidad = {'SVM': list(), 'Rf': list(), 'kNN': list(), 'SVM_Rf': list(), 'SVM_kNN': list(), 'Rf_kNN': list(), 'Prob_media': list(), 'Votos': list()}
    F1_score = {'SVM': list(), 'Rf': list(), 'kNN': list(), 'SVM_Rf': list(), 'SVM_kNN': list(), 'Rf_kNN': list(), 'Prob_media': list(), 'Votos': list()}

    for time, value in cross_val.items():
        
        for fold, val2 in value.items():

            probabilities = {'TL': list(), 'SVM': list(), 'Rf': list(), 'kNN': list()}
            data_train = data.iloc[val2['Train'],:]
            y_train = np.array(y)[val2['Train']]
            data_test = data.iloc[val2['Test'],:]

            # Label true 
            y_true = np.array(y)[val2['Test']]
            probabilities['TL'].extend(y_true)

            # SVM 
            model_svm = SVC(random_state=123, probability=True)
            model_svm.fit(data_train, y_train)
            idx = list(model_svm.classes_).index('E')
            prob_svm = model_svm.predict_proba(data_test)[:, idx].tolist()
            probabilities['SVM'].extend(prob_svm)
                            
            # Random Forest
            model_rf = RandomForestClassifier(random_state=123)
            model_rf.fit(data_train, y_train)
            idx = list(model_rf.classes_).index('E')
            prob_rf = model_rf.predict_proba(data_test)[:, idx].tolist()
            probabilities['Rf'].extend(prob_rf)

            # kNN 
            model_knn = KNeighborsClassifier()
            model_knn.fit(data_train, y_train)
            idx = list(model_knn.classes_).index('E')
            prob_knn = model_knn.predict_proba(data_test)[:, idx].tolist()
            probabilities['kNN'].extend(prob_knn)

            # Sacamos las probabilidades combinadas 
            probabilities['SVM_Rf'] = np.mean([probabilities['SVM'], probabilities['Rf']], axis = 0)
            probabilities['SVM_kNN'] = np.mean([probabilities['SVM'], probabilities['kNN']], axis = 0)
            probabilities['Rf_kNN'] = np.mean([probabilities['Rf'], probabilities['kNN']], axis = 0)
            probabilities['Prob_media'] = np.mean([probabilities['SVM'], probabilities['Rf'], probabilities['kNN']], axis = 0)

            # Rendimiento de los modelos 
            for modelo in accuracy.keys(): 
                if modelo != 'Votos':
                    y_pred = prob_to_clase(probs= probabilities[modelo])
                else:
                    y_pred = modelo_votos(pd.DataFrame(probabilities))
                tn, fp, fn, tp = confusion_matrix(probabilities['TL'], y_pred, labels=['R', 'E']).ravel().tolist()
                accuracy[modelo].append(round((tn+tp)/(tn+tp+fp+fn), 3))
                sensibilidad[modelo].append(round((tp)/(tp+fn), 3))
                especificidad[modelo].append(round(tn/(tn+fp), 3))
                F1_score[modelo].append(round((2*tp)/(2*tp+fp+fn), 3))

    return accuracy, sensibilidad, especificidad, F1_score


def rendimiento_test(data_train, ed_train, data_test, ed_test):
    probabilities = {'sample' : list(data_test.index)}
    prediction = {'sample' : list(data_test.index)}
    performance = pd.DataFrame(np.zeros((4, 8)), 
                           index = ['Acc', 'Sens', 'Spec', 'F1'], 
                           columns = ['SVM', 'Rf', 'kNN', 'Prob_media', 'SVM_Rf', 'SVM_kNN', 'Rf_kNN', 'Votos'])

    # True Label 
    prediction['TL'] = list(ed_test.g_acc)

    # SVM
    model_svm = SVC(random_state=123, probability=True)
    model_svm.fit(data_train, ed_train.g_acc)
    idx = list(model_svm.classes_).index('E')
    prob_svm = model_svm.predict_proba(data_test)[:, idx].tolist()
    probabilities['SVM'] = prob_svm
    prediction['SVM'] = prob_to_clase(prob_svm)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['SVM'], performance = performance, modelo = 'SVM')

    # Random Forest 
    model_rf = RandomForestClassifier(random_state=123)
    model_rf.fit(data_train, ed_train.g_acc)
    idx = list(model_rf.classes_).index('E')
    prob_rf = model_rf.predict_proba(data_test)[:, idx].tolist()
    probabilities['Rf'] = prob_rf
    prediction['Rf'] = prob_to_clase(prob_rf)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['Rf'], performance = performance, modelo = 'Rf')

    # kNN 
    model_knn = KNeighborsClassifier()
    model_knn.fit(data_train, ed_train.g_acc)
    idx = list(model_knn.classes_).index('E')
    prob_knn = model_knn.predict_proba(data_test)[:, idx].tolist()
    probabilities['kNN'] = prob_knn
    prediction['kNN'] = prob_to_clase(prob_knn)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['kNN'], performance = performance, modelo = 'kNN')

    # Combinaciones por probabilidad
    # Probabilidad media
    mean_prob = np.mean(pd.DataFrame(probabilities).loc[:, ['SVM', 'Rf', 'kNN']], axis = 1)
    prediction['Prob_media'] = prob_to_clase(mean_prob)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['Prob_media'], performance = performance, modelo = 'Prob_media')

    # SVM - Rf 
    prob_svm_rf = np.mean(pd.DataFrame(probabilities).loc[:, ['SVM', 'Rf']], axis = 1)
    prediction['SVM_Rf'] = prob_to_clase(prob_svm_rf)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['SVM_Rf'], performance = performance, modelo = 'SVM_Rf')

    # SVM - kNN 
    prob_svm_knn = np.mean(pd.DataFrame(probabilities).loc[:, ['SVM', 'kNN']], axis = 1)
    prediction['SVM_kNN'] = prob_to_clase(prob_svm_knn)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['SVM_kNN'], performance = performance, modelo = 'SVM_kNN')

    # Rf - kNN 
    prob_rf_knn = np.mean(pd.DataFrame(probabilities).loc[:, ['Rf', 'kNN']], axis = 1)
    prediction['Rf_kNN'] = prob_to_clase(prob_rf_knn)
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['Rf_kNN'], performance = performance, modelo = 'Rf_kNN')

    # Votos 
    prediction['Votos'] = modelo_votos(pd.DataFrame(probabilities))
    performance_modelo(y_true = prediction['TL'], y_pred = prediction['Votos'], performance = performance, modelo = 'Votos')

    return pd.DataFrame(performance), pd.DataFrame(probabilities)
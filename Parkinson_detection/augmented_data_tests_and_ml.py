# -*- coding: utf-8 -*-
# augmented_data_tests_and_ML.ipynb

# Automatically generated by Colaboratory.

# Original file is located at
#     https://colab.research.google.com/drive/1EEOCj_KPL8na5NT4WzRCck4FPWVuZEkY
# 

from tensorflow import keras
import tensorflow as tf
from keras.layers import Convolution1D, MaxPooling1D, Dropout, Flatten, Dense,BatchNormalization,UpSampling1D,Dense,Reshape,ZeroPadding1D,AveragePooling1D
from keras.models import Sequential,load_model
import os
import pandas as pd
import numpy as np
import random
from numpy import newaxis
from sklearn import utils
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score,confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu,ranksums
from sklearn.metrics import f1_score,precision_score,recall_score
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings("ignore")

ctgan_custom = './ctgan_custom'
cgan_custom = './ctgan_custom'
ctgan_expanded = './ctgan_expanded'
cgan_expanded = './cgan_expanded'
#directories to load files
cgan_custom_list = os.listdir(cgan_custom)
ctgan_custom_list = os.listdir(ctgan_custom)
cgan_expanded_list = os.listdir(cgan_expanded)
ctgan_expanded_list = os.listdir(ctgan_expanded)

#at first,we check CopulaGAN for 200/2000/20000 custom features(as asked)
#concatenate the dataframes and shuffle to pass into the autoencoder
#try for 200

for file in cgan_custom_list:
  if file == "Healthy_100.csv":
    df_healthy_100 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Healthy_1000.csv":
    df_healthy_1000 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Parcinson_100.csv":
    df_parcinson_100 = pd.read_csv(cgan_custom + '/' + file)
  else:
    df_parcinson_1000 = pd.read_csv(cgan_custom + '/' + file)

df_healthy_100 = df_healthy_100.T[1:]
df_healthy_1000 = df_healthy_1000.T[1:]
df_parcinson_100 = df_parcinson_100.T[1:]
df_parcinson_1000 = df_parcinson_1000.T[1:]
df_healthy_100 = df_healthy_100.T
df_healthy_1000 = df_healthy_1000.T
df_parcinson_100 = df_parcinson_100.T
df_parcinson_1000 = df_parcinson_1000.T

healthy_100 = df_healthy_100.to_numpy()
parcinson_100 = df_parcinson_100.to_numpy()
healthy_1000 = df_healthy_1000.to_numpy()
parcinson_1000 = df_parcinson_1000.to_numpy()

total_200 = np.concatenate((healthy_100,parcinson_100))
total_2000 = np.concatenate((healthy_1000,parcinson_1000))

np.random.shuffle(total_200)
np.random.shuffle(total_2000)

#cgan 1000 samples

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_1000.T[idx], parcinson_1000.T[idx]))
  ranksum.append(ranksums(healthy_1000.T[idx], parcinson_1000.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_cgan_1000 = []
mannwhtn_idx_cgan_1000 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_cgan_1000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_cgan_1000.append(idx)
#print(rank_sum_idx_cgan_1000)
#print(mannwhtn_idx_cgan_1000)

#cgan 100 samples

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_100.T[idx], parcinson_100.T[idx]))
  ranksum.append(ranksums(healthy_100.T[idx], parcinson_100.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_cgan_100 = []
mannwhtn_idx_cgan_100 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_cgan_100.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_cgan_100.append(idx)
#print(rank_sum_idx_cgan_100)
#print(mannwhtn_idx_cgan_100)

for file in ctgan_custom_list:
  if file == "Healthy_100.csv":
    df_healthy_100 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Healthy_1000.csv":
    df_healthy_1000 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Parcinson_100.csv":
    df_parcinson_100 = pd.read_csv(cgan_custom + '/' + file)
  else:
    df_parcinson_1000 = pd.read_csv(cgan_custom + '/' + file)

df_healthy_100 = df_healthy_100.T[1:]
df_healthy_1000 = df_healthy_1000.T[1:]
df_parcinson_100 = df_parcinson_100.T[1:]
df_parcinson_1000 = df_parcinson_1000.T[1:]
df_healthy_100 = df_healthy_100.T
df_healthy_1000 = df_healthy_1000.T
df_parcinson_100 = df_parcinson_100.T
df_parcinson_1000 = df_parcinson_1000.T

healthy_100 = df_healthy_100.to_numpy()
parcinson_100 = df_parcinson_100.to_numpy()
healthy_1000 = df_healthy_1000.to_numpy()
parcinson_1000 = df_parcinson_1000.to_numpy()

total_200 = np.concatenate((healthy_100,parcinson_100))
total_2000 = np.concatenate((healthy_1000,parcinson_1000))

np.random.shuffle(total_200)
np.random.shuffle(total_2000)

#ctgan 1000 samples
mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_1000.T[idx], parcinson_1000.T[idx]))
  ranksum.append(ranksums(healthy_1000.T[idx], parcinson_1000.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_ctgan_1000 = []
mannwhtn_idx_ctgan_1000 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_1000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_1000.append(idx)
#print(rank_sum_idx_ctgan_1000)
#print(mannwhtn_idx_ctgan_1000)

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_100.T[idx], parcinson_100.T[idx]))
  ranksum.append(ranksums(healthy_100.T[idx], parcinson_100.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_ctgan_100 = []
mannwhtn_idx_ctgan_100 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_100.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_100.append(idx)
#print(rank_sum_idx_ctgan_100)
#print(mannwhtn_idx_ctgan_100)


#lower number of samples gives greater number of features passing the tests(for custom datasets)

#same procedure for the expanded datasets

for file in cgan_expanded_list:
  if file == "Healthy_100.csv":
    df_healthy_100 = pd.read_csv(cgan_expanded + '/' + file)
  elif file == "Healthy_1000.csv":
    df_healthy_1000 = pd.read_csv(cgan_expanded + '/' + file)
  elif file == "Parcinson_100.csv":
    df_parcinson_100 = pd.read_csv(cgan_expanded + '/' + file)
  else:
    df_parcinson_1000 = pd.read_csv(cgan_expanded + '/' + file)

df_healthy_100 = df_healthy_100.T[1:]
df_healthy_1000 = df_healthy_1000.T[1:]
df_parcinson_100 = df_parcinson_100.T[1:]
df_parcinson_1000 = df_parcinson_1000.T[1:]
df_healthy_100 = df_healthy_100.T
df_healthy_1000 = df_healthy_1000.T
df_parcinson_100 = df_parcinson_100.T
df_parcinson_1000 = df_parcinson_1000.T

healthy_100 = df_healthy_100.to_numpy()
parcinson_100 = df_parcinson_100.to_numpy()
healthy_1000 = df_healthy_1000.to_numpy()
parcinson_1000 = df_parcinson_1000.to_numpy()

total_200 = np.concatenate((healthy_100,parcinson_100))
total_2000 = np.concatenate((healthy_1000,parcinson_1000))

np.random.shuffle(total_200)
np.random.shuffle(total_2000)

#cgan 1000

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_1000.T[idx], parcinson_1000.T[idx]))
  ranksum.append(ranksums(healthy_1000.T[idx], parcinson_1000.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_cgan_1000 = []
mannwhtn_idx_cgan_1000 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_cgan_1000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_cgan_1000.append(idx)
#print(rank_sum_idx_ctgan_1000)
#print(mannwhtn_idx_ctgan_1000)

#cgan 100

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_100.T[idx], parcinson_100.T[idx]))
  ranksum.append(ranksums(healthy_100.T[idx], parcinson_100.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_cgan_100 = []
mannwhtn_idx_cgan_100 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_cgan_100.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_cgan_100.append(idx)
#print(rank_sum_idx_cgan_100)
#print(mannwhtn_idx_cgan_100)

for file in ctgan_expanded_list:
  if file == "Healthy_100.csv":
    df_healthy_100 = pd.read_csv(ctgan_expanded + '/' + file)
  elif file == "Healthy_1000.csv":
    df_healthy_1000 = pd.read_csv(ctgan_expanded + '/' + file)
  elif file == "Parcinson_100.csv":
    df_parcinson_100 = pd.read_csv(ctgan_expanded + '/' + file)
  else:
    df_parcinson_1000 = pd.read_csv(ctgan_expanded + '/' + file)

df_healthy_100 = df_healthy_100.T[1:]
df_healthy_1000 = df_healthy_1000.T[1:]
df_parcinson_100 = df_parcinson_100.T[1:]
df_parcinson_1000 = df_parcinson_1000.T[1:]
df_healthy_100 = df_healthy_100.T
df_healthy_1000 = df_healthy_1000.T
df_parcinson_100 = df_parcinson_100.T
df_parcinson_1000 = df_parcinson_1000.T

healthy_100 = df_healthy_100.to_numpy()
parcinson_100 = df_parcinson_100.to_numpy()
healthy_1000 = df_healthy_1000.to_numpy()
parcinson_1000 = df_parcinson_1000.to_numpy()

total_200 = np.concatenate((healthy_100,parcinson_100))
total_2000 = np.concatenate((healthy_1000,parcinson_1000))

np.random.shuffle(total_200)
np.random.shuffle(total_2000)

#ctgan 1000

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_1000.T[idx], parcinson_1000.T[idx]))
  ranksum.append(ranksums(healthy_1000.T[idx], parcinson_1000.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_ctgan_1000 = []
mannwhtn_idx_ctgan_1000 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_1000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_1000.append(idx)
#print(rank_sum_idx_ctgan_1000)
#print(mannwhtn_idx_ctgan_1000)

#ctgan 100

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_100.T[idx], parcinson_100.T[idx]))
  ranksum.append(ranksums(healthy_100.T[idx], parcinson_100.T[idx]))
#print(np.asarray(ranksum))
#print(np.asarray(mannwhitney))
rank_sum_idx_ctgan_100 = []
mannwhtn_idx_ctgan_100 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_100.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_100.append(idx)
#print(rank_sum_idx_ctgan_1000)
#print(mannwhtn_idx_ctgan_1000)



# FROM THE ABOVE WE KEEP THE CUSTOM SETS WITH CTGAN METHOD FOR AUGMENTATION
# NOW, WE ARE GOING TO USE RAW FEATURES THAT PASSED THE TESTS AND THEN PCA ON THE SAME FEATURES. SVM TESTED ON BOTH 100 AND 1000 SAMPLES



for file in ctgan_custom_list:
  if file == "Healthy_100.csv":
    df_healthy_100 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Healthy_1000.csv":
    df_healthy_1000 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Parcinson_100.csv":
    df_parcinson_100 = pd.read_csv(cgan_custom + '/' + file)
  elif file == 'Parcinson_1000.csv':
    df_parcinson_1000 = pd.read_csv(cgan_custom + '/' + file)
  elif file == "Healthy_10000.csv":
    df_healthy_10000 = pd.read_csv(cgan_custom + '/' + file)
  else:
    df_parcinson_10000 = pd.read_csv(cgan_custom + '/' + file)

df_healthy_100 = df_healthy_100.T[1:]
df_healthy_1000 = df_healthy_1000.T[1:]
df_parcinson_100 = df_parcinson_100.T[1:]
df_parcinson_1000 = df_parcinson_1000.T[1:]
df_healthy_10000 = df_healthy_10000.T[1:]
df_parcinson_10000 = df_parcinson_10000.T[1:]

df_healthy_100 = df_healthy_100.T
df_healthy_1000 = df_healthy_1000.T
df_parcinson_100 = df_parcinson_100.T
df_parcinson_1000 = df_parcinson_1000.T
df_healthy_10000 = df_healthy_10000.T
df_parcinson_10000 = df_parcinson_10000.T

healthy_100 = df_healthy_100.to_numpy()
parcinson_100 = df_parcinson_100.to_numpy()
healthy_1000 = df_healthy_1000.to_numpy()
parcinson_1000 = df_parcinson_1000.to_numpy()
healthy_10000 = df_healthy_10000.to_numpy()
parcinson_10000 = df_parcinson_10000.to_numpy()

healthy_100_clss = np.zeros(healthy_100.shape[0])
healthy_1000_clss = np.zeros(healthy_1000.shape[0])
healthy_10000_clss = np.zeros(healthy_10000.shape[0])

parcinson_100_clss = np.ones(parcinson_100.shape[0])
parcinson_1000_clss = np.ones(parcinson_1000.shape[0])
parcinson_10000_clss = np.ones(parcinson_10000.shape[0])

total_200 = np.concatenate((healthy_100,parcinson_100))
total_2000 = np.concatenate((healthy_1000,parcinson_1000))
total_20000 = np.concatenate((healthy_10000,parcinson_10000))


total_200_clss = np.concatenate((healthy_100_clss,parcinson_100_clss))
total_2000_clss = np.concatenate((healthy_1000_clss,parcinson_1000_clss))
total_20000_clss = np.concatenate((healthy_10000_clss,parcinson_10000_clss))


total_200,total_200_clss = utils.shuffle(total_200,total_200_clss)
total_2000,total_2000_clss = utils.shuffle(total_2000,total_2000_clss)
total_20000,total_20000_clss = utils.shuffle(total_20000,total_20000_clss)

#ctgan 1000 samples
mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_1000.T[idx], parcinson_1000.T[idx]))
  ranksum.append(ranksums(healthy_1000.T[idx], parcinson_1000.T[idx]))

rank_sum_idx_ctgan_1000 = []
mannwhtn_idx_ctgan_1000 = []

#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)

for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_1000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_1000.append(idx)


#ctgan 10000 samples
mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_10000.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_10000.T[idx], parcinson_10000.T[idx]))
  ranksum.append(ranksums(healthy_10000.T[idx], parcinson_10000.T[idx]))

rank_sum_idx_ctgan_10000 = []
mannwhtn_idx_ctgan_10000 = []
#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_10000.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_10000.append(idx)

#ctgan 100 samples

mannwhitney = []
ranksum = []
alpha = 0.01
for idx in range(healthy_100.shape[1]):
  mannwhitney.append(mannwhitneyu(healthy_100.T[idx], parcinson_100.T[idx]))
  ranksum.append(ranksums(healthy_100.T[idx], parcinson_100.T[idx]))
rank_sum_idx_ctgan_100 = []
mannwhtn_idx_ctgan_100 = []

#check which feature returns true for statistical significance p-val = 1%
#this means, which feature(s) of our dataset have significantly different distribution in our dataset between PD and Healthy(DPR)
for idx in range(np.asarray(ranksum).shape[0]):
  if np.asarray(ranksum)[idx,1]<alpha:
    rank_sum_idx_ctgan_100.append(idx)
  if np.asarray(mannwhitney)[idx,1]<alpha:
    mannwhtn_idx_ctgan_100.append(idx)

def Diff(li1, li2):
    return list(set(li1) - set(li2))

cols = np.zeros(df_healthy_1000.shape[1],dtype = int)
for idx in range(cols.shape[0]):
  cols[idx] = (idx)
print(cols)
cols = cols.tolist()
sig_cols_1000 = Diff(cols,mannwhtn_idx_ctgan_1000)
sig_cols_1000 = (np.asarray(sig_cols_1000))

print(sig_cols_1000)

print(total_2000.shape)
sig_cols_100 = Diff(cols,mannwhtn_idx_ctgan_100)
sig_cols_100 = np.asarray(sig_cols_100)
print(total_200.shape)

# CREATE OUR SUPPORT VECTOR MACHINE(NON-LINEAR)

#grid search for rbf kernel and then apply svm in our data

clf = SVC()

hyperparams = {'gamma':[10,1,0.1,0.01,0.01],
               'C':[100,10,1,0.1,0.01],
               'kernel':['linear'],
               'degree':[2,3,4,5]}

gs = GridSearchCV(clf, param_grid=hyperparams, 
                  scoring="accuracy",
                  n_jobs=-1, cv=10, return_train_score=True)


train_200 = total_200[0:120]
train_200_clss = total_200_clss[0:120]

test_200 = total_200[121:]
test_200_clss = total_200_clss[121:]
gs.fit(train_200, train_200_clss)
print("Mean cross-validated training accuracy score:",
      gs.best_score_)
print("Optimal hyperparameter combination:", gs.best_params_)

preds = gs.best_estimator_.predict(test_200) # Predictions
y_true = test_200_clss # True values


cf_matrix = confusion_matrix(y_true, preds)
print("Test accuracy:", np.round(accuracy_score(y_true, preds), 2))

sns.heatmap(cf_matrix, annot=True, cmap='Blues')
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('True', fontsize=12)

# print(np.sum(train_2000_clss))
# test_2000 = total_2000[1201:]
# test_2000_clss = total_2000_clss[1201:]
# clf.fit(train_2000,train_2000_clss)


#grid search for rbf kernel and then apply svm in our data for 20k samples

clf = SVC()
tree = DecisionTreeClassifier()
nb = GaussianNB()
hyperparams = {#'gamma':[0.1,0.01],
               #'C':[1000,100,10,1,0.1],
               'kernel':['linear'],
               #'degree':[2,3,4,5]
               }

hyper_tree = {'criterion':['entropy','gini'],
              'max_depth':[5,7,9,12],
              'min_samples_split':[2,4,5],
              'min_samples_leaf':[1,2,3,4,5,7,10],
              'splitter':['best'],
              'random_state' : [1,5,10]}
              
hyper_nb = {'var_smoothing':[1e-7,1e-6,1e-5,1e-4,1e-8,1e-9,1e-10,1e-11,1e-12]}

gs = GridSearchCV(nb, param_grid=hyper_nb, 
                  scoring="accuracy",
                  n_jobs=-1, cv=10, return_train_score=True)


train_2000 = total_2000[0:1200]
train_2000_clss = total_2000_clss[0:1200]

test_2000 = total_2000[1201:]
test_2000_clss = total_2000_clss[1201:]
gs.fit(train_2000, train_2000_clss)
print("Mean cross-validated training accuracy score:",
      gs.best_score_)
print("Optimal hyperparameter combination:", gs.best_params_)

preds = gs.best_estimator_.predict(test_2000) # Predictions
y_true = test_2000_clss # True values


cf_matrix = confusion_matrix(y_true, preds)
print("Test accuracy:", np.round(accuracy_score(y_true, preds), 2))

sns.heatmap(cf_matrix, annot=True, cmap='Blues')
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('True', fontsize=12)

# print(np.sum(train_2000_clss))
# test_2000 = total_2000[1201:]
# test_2000_clss = total_2000_clss[1201:]
# clf.fit(train_2000,train_2000_clss)


h = pd.read_csv('./modified_outliers/Healthy.csv')
h = h.T[1:]
h = h.T
h = h.to_numpy()


p = pd.read_csv('./modified_outliers/Parcinson.csv')
p = p.T[1:]
p = p.T
p = p.to_numpy()
TEST = np.concatenate((h,p))
p_c = np.ones(7)
h_c = np.zeros(7)
clss = np.concatenate((h_c,p_c))

preds = gs.best_estimator_.predict(TEST)
cf_matrix = confusion_matrix(clss, preds)

sns.heatmap(cf_matrix, annot=True, cmap='Blues')
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('True', fontsize=12)
from sklearn.metrics import f1_score,precision_score,recall_score
f1 = f1_score(clss,preds)
precision = precision_score(clss,preds)
recall = recall_score(clss,preds)
print("F1 score",f1*100,"%")
print("Precision",precision*100,"%")
print("Recall",recall*100,"%")
print("Accuracy",accuracy_score(clss,preds)*100,"%")

# HERE,USE LESS FEATURES(PASSED OR NOT PASSED OR PCA) AND SAME APPROACH(RBF GRID SEARCH AND LINEAR)

#pca analysis, keep components that their sum is greater than 90%

s = 0
n_components = 2
print(total_20000.shape)
while s <=0.9:
  pca = PCA(n_components=n_components)
  principalComponents = pca.fit_transform(total_20000)
  var = pca.explained_variance_ratio_
  s = np.sum(var)
  n_components+=1
print(s)
print(n_components)
print(principalComponents.shape)

clf = SVC(kernel = 'rbf',gamma = 10,C = 0.01)

hyperparams = {'gamma':[100,10,1,0.1,0.01],
               'C':[100,10,1,0.1,0.01]}

gs = GridSearchCV(clf, param_grid=hyperparams, 
                  scoring="accuracy",
                  n_jobs=-1, cv=10, return_train_score=True)


principalComponents_train = principalComponents[0:12000]
train_2000_clss = total_20000_clss[0:12000]

principalComponents_test= principalComponents[12001:]
test_20000_clss = total_20000_clss[12001:]
gs.fit(principalComponents_train, train_2000_clss)
print("Mean cross-validated training accuracy score:",
      gs.best_score_)
print("Optimal hyperparameter combination:", gs.best_params_)

preds = gs.best_estimator_.predict(principalComponents_test) # Predictions
y_true = test_20000_clss # True values


cf_matrix = confusion_matrix(y_true, preds)
print("Test accuracy:", np.round(accuracy_score(y_true, preds), 2))

sns.heatmap(cf_matrix, annot=True, cmap='Blues')
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('True', fontsize=12)


# print(np.sum(train_2000_clss))
# test_2000 = total_2000[1201:]
# test_2000_clss = total_2000_clss[1201:]
# clf.fit(train_2000,train_2000_clss)


h = pd.read_csv('./dataframes/Healthy.csv')
h = h.T[1:]
h = h.T
h = h.to_numpy()


p = pd.read_csv('./dataframes/Parcinson.csv')
p = p.T[1:]
p = p.T
p = p.to_numpy()
TEST = np.concatenate((h,p))
p_c = np.ones(7)
h_c = np.zeros(7)
clss = np.concatenate((h_c,p_c))


n_test = pca.transform(TEST)

preds = gs.best_estimator_.predict(n_test)
cf_matrix = confusion_matrix(clss, preds)

sns.heatmap(cf_matrix, annot=True, cmap='Blues')
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('True', fontsize=12)
f1 = f1_score(clss,preds)
precision = precision_score(clss,preds)
recall = recall_score(clss,preds)
print("F1 score",f1*100,"%")
print("Precision",precision*100,"%")
print("Recall",recall*100,"%")
print("Accuracy",accuracy_score(clss,preds)*100,"%")

#linear svm

clf = SVC(kernel  = 'linear')
train_200 = principalComponents[0:120]
train_200_clss = total_200_clss[0:120]

test_200 = principalComponents[121:]
test_200_clss = total_200_clss[121:]
clf.fit(train_200,train_200_clss)

preds = clf.predict(test_200)
truths = test_200_clss
cf_matrix = confusion_matrix(truths, preds)
print("Test accuracy:", np.round(accuracy_score(truths, preds), 2))

# sns.heatmap(cf_matrix, annot=True, cmap='Blues')
# plt.xlabel('Predicted', fontsize=12)
# plt.ylabel('True', fontsize=12)


pred = clf.predict(TEST)
cf_matrix = confusion_matrix(clss, pred)
print("Test accuracy:", np.round(accuracy_score(clss, pred), 2))
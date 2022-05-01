# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:48:14 2021

@author: Administor
"""
import numpy as np
import pandas as pd
from keras.models import Model
from keras.layers import Dense, Dropout, Input
from pandas import read_csv

from keras.models import Sequential
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection  import GridSearchCV

seed=7
np.random.seed(seed)

TCGA_sig_feature = pd.read_csv('TCGA_sig_feature.tsv',header='infer',sep="\t")
METABRIC_sig_feature = pd.read_csv('METABRIC_sig_feature.tsv',header='infer',sep="\t")
GSE42568_sig_feature = pd.read_csv('GSE42568_sig_feature.tsv',header='infer',sep="\t")
GSE9893_sig_feature = pd.read_csv('GSE9893_sig_feature.tsv',header='infer',sep="\t")
GSE96058_sig_feature = pd.read_csv('GSE96058_sig_feature.tsv',header='infer',sep="\t")

TCGA_hr = pd.read_csv('sig_HR_TCGA.tsv',header='infer',sep="\t")
METABRIC_hr = pd.read_csv('sig_HR_METABRIC.tsv',header='infer',sep="\t")
GSE42568_hr = pd.read_csv('sig_HR_GSE42568.tsv',header='infer',sep="\t")
GSE9893_hr = pd.read_csv('sig_HR_GSE9893.tsv',header='infer',sep="\t")
GSE96058_hr = pd.read_csv('sig_HR_GSE96058.tsv',header='infer',sep="\t")

TCGA_label = pd.read_csv('TCGA_label.tsv',header='infer',sep="\t")

import math

TCGA_dat = pd.DataFrame()
for i in range(12):
    TCGA_dat[i] = TCGA_sig_feature.iloc[:,i]*math.log(TCGA_hr.iloc[i,0])

METABRIC_dat = pd.DataFrame()
for i in range(12):
    METABRIC_dat[i] = METABRIC_sig_feature.iloc[:,i]*math.log(METABRIC_hr.iloc[i,0])

GSE42568_dat = pd.DataFrame()
for i in range(12):
    GSE42568_dat[i] = GSE42568_sig_feature.iloc[:,i]*math.log(GSE42568_hr.iloc[i,0])

GSE9893_dat = pd.DataFrame()
for i in range(12):
    GSE9893_dat[i] = GSE9893_sig_feature.iloc[:,i]*math.log(GSE9893_hr.iloc[i,0])

GSE96058_dat = pd.DataFrame()
for i in range(12):
    GSE96058_dat[i] = GSE96058_sig_feature.iloc[:,i]*math.log(GSE96058_hr.iloc[i,0])


X = TCGA_dat
Y = TCGA_label

def create_model(optimizer='adam'):
    model = Sequential()
    model.add(Dense(12,activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(8, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(1,activation='sigmoid'))

    model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    return model

model = KerasClassifier(build_fn=create_model)
optimizers = ['adam','SGD']

epochs_1 = list([50,100,150])
batches = list([5,10,20])
param_grid = dict(optimizer=optimizers,epochs=epochs_1,batch_size=batches)

grid = GridSearchCV(estimator=model, param_grid=param_grid)
grid_result = grid.fit(X, Y)
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))

#最佳参数
canshu = grid_result.best_params_
opt = canshu.get('optimizer')
epo = canshu.get('epochs')
bat = canshu.get('batch_size')

##################################
####################################

from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_classification
from sklearn.model_selection import StratifiedKFold


model = KerasClassifier(build_fn=create_model, nb_epoch=epo, batch_size=bat)
#数据标准化，改进算法
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
steps = []
steps.append(('standardize', StandardScaler()))
steps.append(('mlp', model))
pipeline = Pipeline(steps)

kfold = KFold(n_splits=10, shuffle=True, random_state=seed)
cv1= StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
results = cross_val_score(pipeline, X, Y, cv=kfold)
results = cross_val_score(pipeline, X, Y, cv=cv1)
print('Standardize: %.2f (%.2f) MSE' % (results.mean(), results.std()))
print(results.mean())
result1 = pd.DataFrame(results)
result_mean = results.mean()
result1.iloc[9,0] = result_mean
result1.index = [1,2,3,4,5,6,7,8,9,'mean']
result1.columns = ['cross_val_score']
result1.to_csv('cross_val_score.tsv',sep="\t")


# Hisotry列表
print(results.history.keys())

# accuracy的历史
from matplotlib import pyplot as plt
plt.plot(results.history['acc'])
plt.plot(results.history['val_acc'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()

# loss的历史
plt.plot(results.history['loss'])
plt.plot(results.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'validation'], loc='upper left')
plt.show()
######################################################
######################ROC
import pandas as pd
from sklearn.metrics import precision_recall_curve, auc, roc_curve
import math
import matplotlib.pyplot as plt

X = np.array(X)
Y = TCGA_label
Y = np.array(Y)
Y[Y == 1] = 0
Y[Y == 2] = 1

axes = plt.subplots(1,1,figsize=(6, 6),dpi=300)

cv1= StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
predictor = KerasClassifier(build_fn=create_model, nb_epoch=epo, batch_size=bat)

y_real = []
y_proba = []
for i, (train_index, test_index) in enumerate(cv1.split(X,Y)):
    Xtrain, Xtest = X[train_index], X[test_index]
    Ytrain, Ytest = Y[train_index], Y[test_index]
    predictor.fit(Xtrain, Ytrain)
    pred_proba = predictor.predict_proba(Xtest)
    
    #y_pred_keras = predictor.predict(Xtest).ravel()
    #fpr_keras, tpr_keras, thresholds_keras = roc_curve(Ytest, y_pred_keras)
    #precision, recall, _ = precision_recall_curve(Ytest, pred_proba[:,1])
    #lab = 'Fold %d AUC=%.4f' % (i+1, auc(recall, precision))
    #axes[1].step(recall, precision, label=lab)
    y_real.append(Ytest)
    y_proba.append(pred_proba[:,1])

y_real = np.concatenate(y_real)
y_proba = np.concatenate(y_proba)
precision, recall, _ = precision_recall_curve(y_real, y_proba)
lab = 'Overall AUC=%.4f' % (auc(recall, precision))


axes[1].step((1-precision),recall,label=lab, lw=2, color='red')
axes[1].set_title('ROC curve',fontsize =12)
axes[1].set_xlabel('1-Precision',fontsize =12)
axes[1].set_ylabel('Recall',fontsize =12)
axes[1].legend(loc='lower right', fontsize=12)
#axes[1].legend(loc='upper left', fontsize=12)
plt.show()

axes[0].savefig('result.pdf')
###########################################
# make predictions
from sklearn.preprocessing import StandardScaler
estimator = KerasClassifier(build_fn=create_model, epochs=epo, batch_size=bat)
estimator.fit(X, Y)
# make predictions
METABRIC_scale= StandardScaler().fit_transform(METABRIC_dat)
pred = estimator.predict(METABRIC_scale)
pred_METABRIC = pd.DataFrame(pred)
pred_METABRIC.to_csv('METABRIC_predict_lable.tsv',sep="\t",header=False)

GSE42568_scale= StandardScaler().fit_transform(GSE42568_dat)
pred_GSE42568 = estimator.predict(GSE42568_scale)
pred_GSE42568 = pd.DataFrame(pred_GSE42568)
pred_GSE42568.to_csv('GSE42568_predict_lable.tsv',sep="\t",header=False)

GSE9893_scale= StandardScaler().fit_transform(GSE9893_dat)
new_GSE9893 = StandardScaler().fit_transform(GSE9893_scale)
pred_GSE9893 = estimator.predict(new_GSE9893)
pred_GSE9893 = pd.DataFrame(pred_GSE9893)
pred_GSE9893.to_csv('GSE9893_predict_lable.tsv',sep="\t",header=False)

GSE96058_scale= StandardScaler().fit_transform(GSE96058_dat)
pred_GSE96058 = estimator.predict(GSE96058_scale)
pred_GSE96058 = pd.DataFrame(pred_GSE96058)
pred_GSE96058.to_csv('GSE96058_predict_lable.tsv',sep="\t",header=False)

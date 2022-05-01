# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 09:28:59 2021

@author: Administor
"""

import numpy as np
import pandas as pd
from keras.models import Model
from keras.layers import Dense, Dropout, Input
from keras import regularizers
#设置随机种子
seed = 7
np.random.seed(seed)

train_data = pd.read_csv('BRCA_TF_train.tsv',sep="\t")
test_data1 = pd.read_csv('BRCA_TF_test1.tsv',sep="\t")
test_data2 = pd.read_csv('BRCA_TF_test2.tsv',sep="\t")
test_data3 = pd.read_csv('BRCA_TF_test3.tsv',sep="\t")
test_data4 = pd.read_csv('BRCA_TF_test4.tsv',sep="\t")

x_train = train_data.T
x_train = x_train.astype('float64')

x_test1 = test_data1.T
x_test1 = x_test1.astype('float64')

x_test2 = test_data2.T
x_test2 = x_test2.astype('float64')

x_test3 = test_data3.T
x_test3 = x_test3.astype('float64')

x_test4 = test_data4.T
x_test4 = x_test4.astype('float64')

#将数据标准到0-1之间
x_train = x_train.apply(lambda x:(x - np.min(x))/(np.max(x) - np.min(x)),axis=0)
x_test1 = x_test1.apply(lambda x:(x - np.min(x))/(np.max(x) - np.min(x)),axis=0)
x_test2 = x_test2.apply(lambda x:(x - np.min(x))/(np.max(x) - np.min(x)),axis=0)
x_test3 = x_test3.apply(lambda x:(x - np.min(x))/(np.max(x) - np.min(x)),axis=0)
x_test4 = x_test4.apply(lambda x:(x - np.min(x))/(np.max(x) - np.min(x)),axis=0)


bottleneck_size=100
input_size = 1001
hidden_size =500
#code_size = 1000

input_data = Input(shape=(input_size,))
## encoded部分
encoded = Dense(hidden_size,activation="relu",
                kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                bias_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                activity_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-6),
                name = "enco1")(input_data)
Dropout(0.5)
#encoded = Dense(code_size,activation="tanh",name = "enco2")(encoded)
#Dropout(0.2)
encoded = Dense(bottleneck_size,activation="tanh",
                kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                bias_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                activity_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-6),
                name = "enco2")(encoded)
Dropout(0.5)
## decoded部分

decoded = Dense(hidden_size,activation="relu",
                kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                bias_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                activity_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-6),
                name = "deco1")(encoded)
Dropout(0.5)
decoded = Dense(input_size,activation="tanh",
                kernel_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                bias_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-5),
                activity_regularizer=regularizers.l1_l2(l1=1e-5, l2=1e-6),
                name = "deco2")(decoded)
Dropout(0.5)

## 连接为自编码模型
autoencoder = Model(input_data, decoded)
autoencoder.summary()
autoencoder.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
## 模型训练
autoencoder_fit = autoencoder.fit(x_train, x_train,epochs=100,
                                  shuffle=True,validation_data=(x_test1,x_test1))

#autoencoder_data = autoencoder.predict(x_train)

## 先定义encoder模型作为输出
encoder = Model(input_data, encoded)
encoder.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
## 获取训练集的encoded
encoded_data = encoder.predict(x_train)
colna = ["feature"+str(ii+1) for ii in range(bottleneck_size)]
train_encoded_feture = pd.DataFrame(data = encoded_data,columns=colna)
train_encoded_feture.index = x_train.index

test_encoded_data1 = encoder.predict(x_test1)
test_encoded_feture1 = pd.DataFrame(data = test_encoded_data1,columns=colna)
test_encoded_feture1.index = x_test1.index

test_encoded_data2 = encoder.predict(x_test2)
test_encoded_feture2 = pd.DataFrame(data = test_encoded_data2,columns=colna)
test_encoded_feture2.index = x_test2.index

test_encoded_data3 = encoder.predict(x_test3)
test_encoded_feture3 = pd.DataFrame(data = test_encoded_data3,columns=colna)
test_encoded_feture3.index = x_test3.index

test_encoded_data4 = encoder.predict(x_test4)
test_encoded_feture4 = pd.DataFrame(data = test_encoded_data4,columns=colna)
test_encoded_feture4.index = x_test4.index


train_encoded_feture.to_csv('BRCAtrain_TCGA_encoded_feature.tsv',sep="\t")
test_encoded_feture1.to_csv('BRCAtest1_GSE42568_encoded_feature.tsv',sep="\t")
test_encoded_feture2.to_csv('BRCAtest2_GSE9893_encoded_feature.tsv',sep="\t")
test_encoded_feture3.to_csv('BRCAtest3_METABRIC_encoded_feature.tsv',sep="\t")
test_encoded_feture4.to_csv('BRCAtest4_GSE96058_encoded_feature.tsv',sep="\t")
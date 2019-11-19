# Author: Uygar Sumbul, Olga Gliko, Rohan Gala
# Allen Institute

import numpy as np
import keras
import scipy as sp
import scipy.io as sio
from scipy.stats import norm
from keras.layers import Input, Dense, Lambda, Layer, Dropout, BatchNormalization
from keras.models import Model
from keras import backend as K
from keras import metrics
from keras.objectives import binary_crossentropy
from keras.callbacks import LearningRateScheduler
from keras.losses import mean_squared_error, mean_absolute_error
from keras.regularizers import l2
from keras.constraints import unit_norm

import tensorflow as tf
import sys
import os

# load data from .mat file containing expression of 6083 highly expressed genes 
# and 47 neuropeptide genes in 22,439 neurons along with sample id
root_folder              = "/nas5/peptides"
data                     = sio.loadmat(os.path.join(root_folder, 'mouse_V1_ALM_20180520_6_5byExpression_and_NP18andGPCR29.mat'))
pep                      = data['pep']
sample_id                = data['sample_id']
thisRun                  = int(sys.argv[2])
foldCount                = 13
foldSize                 = pep.shape[0] / foldCount
heldOutInd               = (np.arange(thisRun*foldSize, (thisRun+1)*foldSize)).astype('int')
trainingInd              = (np.setdiff1d(np.arange(pep.shape[0]), heldOutInd)).astype('int')
heldOutPep               = pep[heldOutInd,  :]
pep                      = pep[trainingInd, :]
original_dim_pep         = pep.shape[1]
intermediate_dim1        = 100
intermediate_dim2        = 50%64
bottleneck_dim           = int(sys.argv[1])
n_epoch1                 = 10000
bat_size                 = 794
dropoutRate1             = 0.8
dropoutRate2             = 0.0
full_train_size          = trainingInd.size
bb                       = int(sys.argv[1])
ff                       = trainingInd.size

# load a latent space representation of single autoencoder from .mat file
data = sio.loadmat(os.path.join(root_folder, 'singleAE_6_5byExpression_dim5_run0_iter50K_0.8Dropout_intermediate100_BN_bat956.mat'))
z1 = data['e1']
val_z1 = data['et1']

y               = Input(shape=(original_dim_pep,),               name='y')
hidden2         = Dropout(dropoutRate2, name='drop2')(y)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense10')(hidden2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense11')(hidden2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense12')(hidden2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense13')(hidden2)
hidden2         = Dense(bottleneck_dim,    activation='linear', name='dense14')(hidden2)
z2              = BatchNormalization(name='z2', center=False, scale=False,epsilon=1e-10)(hidden2)
hidden2         = Dense(intermediate_dim2, activation='linear', name='dense15')(z2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense16')(hidden2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense17')(hidden2)
hidden2         = Dense(intermediate_dim2, activation='relu', name='dense18')(hidden2)
xd2             = Dense(original_dim_pep,  activation='relu',    name='xd2')(hidden2)

caeDual         = Model(inputs=[y], outputs=[xd2, z2])

def cae1_loss(y_true, y_pred):
  return mean_squared_error(y_true, y_pred)
def loss_latentdim(y_true, y_pred):
  zz1   = y_true - tf.reduce_mean(y_true, axis=0)
  zz2   = y_pred - tf.reduce_mean(y_pred, axis=0)
  s1    = tf.svd(zz1, compute_uv=False)
  mins1 = tf.reduce_min(tf.square(s1))
  s2    = tf.svd(zz2, compute_uv=False)
  mins2 = tf.reduce_min(tf.square(s2))
  denom = tf.minimum(mins1,mins2)
  C = mean_squared_error(y_true, y_pred)/denom
  return C
                                                                                                                                            
caeDual.compile(optimizer='adam', loss={'xd2': cae1_loss, 'z2': loss_latentdim}, loss_weights={'xd2': 1., 'z2': 100.})
history         = caeDual.fit({'y':pep}, {'xd2' : pep, 'z2' : z1}, batch_size=bat_size, epochs=n_epoch1, 
                               validation_data=({'y':heldOutPep}, {'xd2':heldOutPep, 'z2':val_z1}))

result          = caeDual.predict(pep)
e2              = result[1]
result          = caeDual.predict(heldOutPep)
et2             = result[1]
val_xd2_loss    = history.history['val_xd2_loss']
xd2_loss        = history.history['xd2_loss']
val_xd2_loss    = val_xd2_loss[::10]
xd2_loss        = xd2_loss[::10]
fileName        = os.path.join(root_folder, 'dualAE_inputZ1_6_5byExpression_and_NP18GPCR29_dim' + sys.argv[1] 
                               + '_run' + sys.argv[2] + '_iter10K_loss1_100_0.0Dropout_intermediate_50_bat794_neuronsOnly.mat')
sio.savemat(fileName, {'e2':e2, 'et2':et2, 'sample_id':sample_id, 'val_xd2_loss':val_xd2_loss, 'xd2_loss':xd2_loss})

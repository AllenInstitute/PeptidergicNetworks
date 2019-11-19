# Authors: Uygar Sumbul, Olga Gliko, Rohan Gala
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
data                     = sio.loadmat(os.path.join(root_folder,'mouse_V1_ALM_20180520_6_5byExpression_and_NP18andGPCR29.mat'))
logOnePlusGeneExpression = data['logOnePlusGeneExpression']
sample_id                = data['sample_id']
thisRun                  = int(sys.argv[2])
foldCount                = 13
foldSize                 = logOnePlusGeneExpression.shape[0] / foldCount
heldOutInd               = (np.arange(thisRun*foldSize, (thisRun+1)*foldSize)).astype('int')
trainingInd              = (np.setdiff1d(np.arange(logOnePlusGeneExpression.shape[0]), heldOutInd)).astype('int')
heldOut                  = logOnePlusGeneExpression[heldOutInd,  :]
logOnePlusGeneExpression = logOnePlusGeneExpression[trainingInd, :]
lbl                      = (logOnePlusGeneExpression>0).astype(np.float32)
original_dim             = logOnePlusGeneExpression.shape[1]
intermediate_dim1        = 100
bottleneck_dim           = int(sys.argv[1])
n_epoch1                 = 50000
bat_size                 = 956
dropoutRate1             = 0.8
full_train_size          = trainingInd.size
bb                       = int(sys.argv[1])
ff                       = trainingInd.size

x               = Input(shape=(original_dim,),                   name='x')
hidden1         = Dropout(dropoutRate1, name='drop1')(x)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense1')(hidden1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense2')(hidden1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense3')(hidden1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense4')(hidden1)
hidden1         = Dense(bottleneck_dim,    activation='linear', name='dense5')(hidden1)
z1              = BatchNormalization(name='z1', center=False, scale=False,epsilon=1e-10)(hidden1)
hidden1         = Dense(intermediate_dim1, activation='linear', name='dense6')(z1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense7')(hidden1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense8')(hidden1)
hidden1         = Dense(intermediate_dim1, activation='relu', name='dense9')(hidden1)
xd1             = Dense(original_dim,      activation='relu', name='xd1')(hidden1)

caeDual         = Model(inputs=x, outputs=[xd1, z1])
def cae1_loss(y_true, y_pred):
  return mean_squared_error(y_true, y_pred)

caeDual.compile(optimizer='adam', loss={'xd1': cae1_loss, 'z1': cae1_loss}, loss_weights={'xd1': 1., 'z1': 0.})
history         = caeDual.fit({'x' : logOnePlusGeneExpression}, {'xd1' : logOnePlusGeneExpression, 
                              'z1' : np.zeros((ff, bb))}, batch_size=bat_size, epochs=n_epoch1, 
                              validation_data=(heldOut, {'xd1':heldOut, 'z1':np.zeros((heldOut.shape[0], bb))}))

result          = caeDual.predict(logOnePlusGeneExpression)
e1              = result[1]
d1              = result[0]
result          = caeDual.predict(heldOut)
et1             = result[1]
dt1             = result[0]
val_xd1_loss    = history.history['val_xd1_loss']
xd1_loss        = history.history['xd1_loss']
val_xd1_loss    = val_xd1_loss[::10]
xd1_loss        = xd1_loss[::10]
fileName        = os.path.join(root_folder, 'singleAE_6_5byExpression_dim' + sys.argv[1] + '_run' + sys.argv[2] 
                               + '_iter50K_0.8Dropout_intermediate100_BN_bat956.mat')
sio.savemat(fileName, {'e1':e1, 'd1':d1, 'et1':et1, 'dt1':dt1, 'sample_id':sample_id, 'val_xd1_loss':val_xd1_loss, 
            'xd1_loss':xd1_loss})
modelweights_filename = os.path.join(root_folder, 'singleAE_model0004_weights.h5')
caeDual.save_weights(modelweights_filename)
fileName              = os.path.join(root_folder, 'kerasModel_singleAE_6_5byExpression_dim' + sys.argv[1] + '_run' 
                                     + sys.argv[2] + '_iter50K_0.8Dropout_intermediate100_BN_bat956.h5')
caeDual.save(fileName)




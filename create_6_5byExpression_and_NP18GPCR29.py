# Authors: Uygar Sumbul, Olga Gliko, Rohan Gala
# Allen Institute

import os
import feather
import numpy as np
import scipy.io as sio

# load data from .feather files containing expression of 21,931 genes in 22,439 neurons 
# along with sample id, class label, etc.
root_folder = "/nas5/peptides"
path               = os.path.join(root_folder,'mouse_V1_ALM_20180520_data.feather')
df                 = feather.read_dataframe(path)
path               = os.path.join(root_folder,'mouse_V1_ALM_20180520_anno.feather')
dfanno             = feather.read_dataframe(path)
df                 = df[dfanno['class_label'].isin(['GABAergic', 'Glutamatergic'])]
dfanno             = dfanno[dfanno['class_label'].isin(['GABAergic', 'Glutamatergic'])]
lOPGE              = np.log(1 + np.array(df.values[:,1:], dtype=np.float32))
sample_id          = np.array(df.loc[:,'sample_id'], dtype=np.object)
gene_id            = df.columns[1:]
# specify 47 neuropeptide genes
gene_id_pep        = ['Vip', 'Npy', 'Sst', 'Penk', 'Tac2', 'Cck', 'Crh', 'Tac1', 'Pdyn', 'Cort', 'Pthlh', 'Pnoc',
                      'Adcyap1', 'Trh', 'Grp', 'Nmb', 'Nts', 'Rln1', 'Vipr1', 'Vipr2', 'Npy1r', 'Npy2r', 'Npy5r',
                      'Sstr1', 'Sstr2', 'Sstr3', 'Sstr4', 'Oprd1', 'Oprm1', 'Tacr3', 'Cckbr', 'Crhr1', 'Crhr2', 
                      'Tacr1', 'Oprk1', 'Pth1r', 'Oprl1', 'Adcyap1r1', 'Trhr', 'Trhr2', 'Grpr', 'Nmbr', 'Ntsr1',
                      'Ntsr2', 'Rxfp1', 'Rxfp2', 'Rxfp3']
pep                = np.log(1 + np.array(df[gene_id_pep], dtype=np.float32))
# select highly expressed genes with log(Expression + 1)>6.5
richExpression     = np.sum((lOPGE>6.5).astype(np.int), axis=0)
richExpressionLocs = richExpression>0
lOPGE              = lOPGE[:, richExpressionLocs]
gene_id            = np.array(gene_id[richExpressionLocs],              dtype=np.object)
cre                = np.array(dfanno.loc[:,'cre_label'],                dtype=np.object)
cluster            = np.array(dfanno.loc[:,'cluster_label'],            dtype=np.object)
core               = np.array(dfanno.loc[:,'core_int_label'],           dtype=np.object)
genotype           = np.array(dfanno.loc[:,'genotype_label'],           dtype=np.object)
layer              = np.array(dfanno.loc[:,'layer_label'],              dtype=np.object)
c1                 = np.array(dfanno.loc[:,'primary_cluster_label'],    dtype=np.object)
c2                 = np.array(dfanno.loc[:,'secondary_cluster_label'],  dtype=np.object)
creID              = np.array(dfanno.loc[:,'cre_id'],                   dtype=np.object)
clusterID          = np.array(dfanno.loc[:,'cluster_id'],               dtype=np.object)
coreID             = np.array(dfanno.loc[:,'core_int_id'],              dtype=np.object)
genotypeID         = np.array(dfanno.loc[:,'genotype_id'],              dtype=np.object)
layerID            = np.array(dfanno.loc[:,'layer_id'],                 dtype=np.object)
c1ID               = np.array(dfanno.loc[:,'primary_cluster_id'],       dtype=np.object)
c2ID               = np.array(dfanno.loc[:,'secondary_cluster_id'],     dtype=np.object)
sample_id_anno     = np.array(dfanno.loc[:,'sample_id'],                dtype=np.object)
region             = np.array(dfanno.loc[:,'region_label'],             dtype=np.object)
shuffledInd        = np.arange(lOPGE.shape[0])
np.random.seed(0)
np.random.shuffle(shuffledInd)
lOPGE              = lOPGE[shuffledInd, :]
pep                = pep[shuffledInd, :]
sample_id          = sample_id[shuffledInd]
cre                = cre[shuffledInd]
cluster            = cluster[shuffledInd]
core               = core[shuffledInd]
genotype           = genotype[shuffledInd]
layer              = layer[shuffledInd]
c1                 = c1[shuffledInd]
c2                 = c2[shuffledInd]
creID              = creID[shuffledInd]
clusterID          = clusterID[shuffledInd]
coreID             = coreID[shuffledInd]
genotypeID         = genotypeID[shuffledInd]
layerID            = layerID[shuffledInd]
c1ID               = c1ID[shuffledInd]
c2ID               = c2ID[shuffledInd]
sample_id_anno     = sample_id_anno[shuffledInd]
region             = region[shuffledInd]
path               = os.path.join(root_folder,'mouse_V1_ALM_20180520_6_5byExpression_and_NP18andGPCR29_region.mat')
sio.savemat(path,{'logOnePlusGeneExpression':lOPGE, 'sample_id':sample_id, 'gene_id':gene_id, 'cre':cre,
            'cluster':cluster, 'core':core, 'genotype':genotype, 'layer':layer, 'c1':c1,'c2':c2, 'creID':creID,
            'clusterID':clusterID, 'coreID':coreID, 'genotypeID':genotypeID, 'layerID':layerID,
            'sample_id_anno':sample_id_anno, 'region':region, 'pep':pep, 'gene_id_pep':gene_id_pep,
            'shuffledInd':shuffledInd})

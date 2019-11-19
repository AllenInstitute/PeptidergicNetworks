import numpy as np 
import scipy
import scipy.io as sio
import sklearn
from scipy.stats import multivariate_normal as mvn
import os

log_file = "/nas5/pepetides/dimensionality_testing.json"

f = open(log_file,"w+")
supervised_params = {
    "matfile": matfile,
    "ground_truth": ground_truth,
    "num_clusters": clusterID.shape[0],
    "foldcount":foldCount,
    "supervised_ari": ari
}

unsupervised_params = {
    "matfile": matfile,
    "num_clusters": clusterID.shape[0],
    "unsupervised_ari": ari
    "GMM_params": {
        "repeats":repeats,
        "replicates":,
        "regularization":,
    }
}
########################################
####Supervised GMM fit to encoding####
########################################
matfile = "singleAE_BN_15092018_6byExpression_dataOct31_annoNov28_dim2_run0_iter20K_0.8Dropout_intermediate100_BN_bat956repeat.mat"
ground_truth = "/nas5/peptides/aibs_mouse_facsseq_v1_alm_20170913_6byExpression_and_NP19andGPCR24_dataOct31_annoNov28.mat"
mat_contents = sio.loadmat(matfile)
gt_contents = sio.loadmat(ground_truth)
clusterID = np.concatenate(gt_contents["clusterID"])
e1 =mat_contents["e1"]
thisRun = 0
foldCount = 13
foldSize = gt_contents["logOnePlusGeneExpression"].shape[0] / foldCount
all_ind = np.arange(gt_contents["logOnePlusGeneExpression"].shape[0])
val_ind = np.arange(foldSize)
train_ind  = np.asarray(list(set(all_ind) - set(val_ind)))
cluster_num = int(max(clusterID))
clusterLik = np.zeros((mat_contents["et1"].size , cluster_num))
clusterMeans = np.empty(clusterID.size,2)
clusterCovs = np.empty((clusterID.size, train_ind.size))
for cluster in clusterID: 
    clusterSamples = np.squeeze(np.array([e1[c-val_ind.shape[0]] for c in train_ind if clusterID[c] == cluster]))
    clusterMeans[int(cluster),:] = np.mean(clusterSamples, axis =0)
    if clusterSamples.shape[0]<50:
        clusterCovs[int(cluster)] = np.diagonal(np.diagonal(np.cov(np.squeeze(clusterSamples))))
    else:
        clusterCovs[int(cluster),:] = np.cov(clusterSamples,rowvar=False)
    clusterLik[:, int(cluster)] = mvn.pdf(mat_contents["et1"], clusterMeans[int(cluster)], clusterCovs[int(cluster)])
valClusterID = clusterID[:val_ind.shape[0]]
clusterAssign = np.argmax(clusterLik,axis=1)
np.count_nonzero(clusterAssign)
#disp(nnz(clusterAss(:)==valClusterID(:))/numel(clusterAss)); % 0.8894
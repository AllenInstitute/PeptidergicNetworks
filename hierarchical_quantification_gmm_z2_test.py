# Raw hierarchical tree file: /data/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData
# tree_20180520.csv was obtained using `dend_functions.R` and `dend_parents.R` functions (Ref. Rohan/Zizhen)

# To visually inspect the original tree:
# http://molgen-shiny.corp.alleninstitute.org/heatmap_2017/?db=/data/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520

import matplotlib
matplotlib.use('Agg')
from scipy.stats import multivariate_normal as mvn
import numpy as np
import scipy.io as sio
import feather
import pandas as pd
import datetime
import os
import csv
import matplotlib.pyplot as plt

def merge(y,child,parent,clusterID,lookup,merged_nodes):
    minind = np.nanargmin(y)
    this_parent = child[minind]
    c_ind = np.where(parent==this_parent)[0]
    lookup[this_parent]=max(lookup.values())+1
    for i in c_ind:
        l = child[i]
        #print(i,': ',child[i],'-->',this_parent)
        if lookup.has_key(l):
            old_clusterID = lookup[l]
            old_ids = np.where(clusterID==old_clusterID)[0]
            
            #print(i,': ',child[i],'-->',this_parent)
            child[i]=this_parent
            y[minind]=np.nan
            
            #print(old_clusterID,'-->', lookup[this_parent])
            for oi in old_ids:
                clusterID[oi]=lookup[this_parent]
        else:
            #removing glia from dictionary and setting height to none
            lookup.pop(this_parent, None)
            child[i]=this_parent
            y[minind]=np.nan
    #keeping track of merged nodes
    merged_nodes.append(this_parent)
    return merged_nodes


def calc_rindex(predict_history,heights,children,neuron_h):
    res_indices =np.ones(predict_history.shape[0])
    for cell_i in range(predict_history.shape[0]):
        correct = np.where(predict_history[cell_i,:]+1==clusterID_history[cell_i,:])[0]
        lowest = np.min(correct)
        for l in lookup.keys():
            if lookup[l]== clusterID_history[cell_i,lowest]:
                #print('label found: ',l)
                ind = np.where(children == l)
        res_indices[cell_i] = heights[ind]
    res = np.mean(np.array([1-i for i in res_indices]))
    #res_norm = np.mean(np.array([1-(i/neuron_h) if i <= neuron_h else 0 for i in res_indices]))
    return res,res_indices

start = datetime.datetime.now()
htree = pd.read_csv('/nas5/peptides/tree_20180520.csv')
htree = htree[['x','y','leaf','label','parent','col']]
htree['leaf'] = htree['leaf'].values==True
htree = htree.sort_values(by=['y','x'], axis=0, ascending=[True,True]).copy(deep=True)
htree = htree.reset_index(drop=True).copy(deep=True)


ground_truth = "/nas5/peptides/mouse_V1_ALM_20180520_6_5byExpression_and_NP18andGPCR29.mat"
gt_contents = sio.loadmat(ground_truth)
clusterID = np.squeeze(np.concatenate(np.concatenate(gt_contents["clusterID"])))
clusters = np.squeeze(np.concatenate(np.concatenate(gt_contents["cluster"])))

#Copy of fields used for merging:
y = htree['y'].values.copy()
y[htree['leaf'].values]=np.nan
child = htree['label'].values.copy()
parent = htree['parent'].values.copy()
leaf = htree['leaf'].values.copy()
#for res_score
children = htree['label'].values.copy()
heights = htree['y'].values.copy() #for res_score
clusterID = np.squeeze(np.concatenate(np.concatenate(gt_contents["clusterID"])))
clusters = np.squeeze(np.concatenate(np.concatenate(gt_contents["cluster"])))
matfile = "/nas5/peptides/dualAE_inputZ1_6_5byExpression_and_NP18GPCR29_dim5_run0_iter10K_loss1_100_0.0Dropout_intermediate_50_bat794_neuronsOnly.mat"
mat_contents = sio.loadmat(matfile)
split = matfile.split('_')
np_set = "dataOct31"
for d in split:
    if "dim" in d:
        dims = d.split("m")
        e = dims[1] 
    elif "NP" in d:
        np_set = d
#start = datetime.datetime.now()

#dict to keep track of cell-type labels and ground truth clusterID
lookup = {}
for c, id in zip(clusters,clusterID):
    lookup[str(c)]= id

e2 = np.array(mat_contents["e2"])
et2 = np.array(mat_contents["et2"])
e2_all = np.concatenate((et2,e2))
thisRun = 0
foldCount = 13
foldSize = gt_contents["logOnePlusGeneExpression"].shape[0] / foldCount
all_ind = np.arange(gt_contents["logOnePlusGeneExpression"].shape[0])
val_ind = np.arange(foldSize+1)
train_ind  = np.asarray(list(set(all_ind) - set(val_ind)))
cuts = np.where(leaf ==False )[0]
predict_history = np.zeros((val_ind.shape[0],cuts.shape[0]))
clusterID_history = np.zeros((val_ind.shape[0],cuts.shape[0]))

cut=0
merged_nodes = [] 

while cut < cuts.shape[0]:
    cluster_num = int(max(clusterID))
    clusterLik = np.zeros((mat_contents["et2"].shape[0], cluster_num))
    clusterMeans = np.empty((cluster_num,int(e)))
    clusterCovs = [None]*(cluster_num) 

    for cluster in range(cluster_num):
        clusterSamples = e2[clusterID[train_ind]==cluster+1,:]
        if clusterSamples.shape[0]>0:
            means = np.mean(clusterSamples, axis =0)
            clusterMeans[int(cluster),:] = means
            if clusterSamples.shape[0]<50:
                clusterCovs[int(cluster)] = np.diagonal(np.cov(clusterSamples,rowvar=False))
            else:
                clusterCovs[int(cluster)] = np.cov(clusterSamples,rowvar=False)

            clusterLik[:, int(cluster)] = mvn.pdf(et2,clusterMeans[int(cluster),:], 
                                                         clusterCovs[int(cluster)])
    
    valClusterID = clusterID[:val_ind.shape[0]]
    clusterID_history[:,cut] = valClusterID
    clusterAssign = clusterLik.argmax(1)
    predict_history[:,cut] = clusterAssign


    #predict_history = validate(clusterID,clusterLik,predict_history,cut)
    merged_nodes = merge(y,child,parent,clusterID,lookup,merged_nodes)
    cut = cut +1

#grabbing height of n4, first neural node
neuron = np.where(children == 'n3')
neuron_h = heights[neuron]
norm = 1 - neuron_h
#calculating avg resIndex 
res, res_indices = calc_rindex(predict_history,heights,children,neuron_h)
res_norm = float((res-norm)/(1-norm)) #normalizing over n3 (1-0.46)
end = datetime.datetime.now()
sup_time = end - start
print(res_norm, res)

sup_log_file = "/nas5/peptides/retesting.csv"
if not os.path.exists(sup_log_file):
    open(sup_log_file,'a').close()

tree_params = {
    "run": 10,
    "nodes":cuts.shape[0],
    "matfile": matfile,
    "ground_truth": ground_truth,
    "heights": None,
    "norm_node": 'n3',
    "res_index": res,
    "norm_res_index": res_norm,
    "time_to_run": str(sup_time),
    "z": "z2"
}

with open(sup_log_file, 'a') as csvfile:
    headers = ["run","norm_node", "z", "res_index","nodes","heights","time_to_run", 
                            "norm_res_index",
                            "matfile","ground_truth"]
    writer = csv.DictWriter(csvfile,fieldnames = headers)
    writer.writerow(tree_params)


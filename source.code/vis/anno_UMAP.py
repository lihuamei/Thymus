import os, sys, scipy
import time, pdb
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from bbknn import bbknn

def regress_batch_v2(adata,batch_key,confounder_key):
    '''batch regression tool
    batch_key=list of observation categories to be regressed out
    confounder_key=list of observation categories to be kept
    returns ndata with corrected X'''

    from sklearn.linear_model import Ridge
    dummy = pd.get_dummies(adata.obs[batch_key+confounder_key],drop_first=False)
    X_exp = adata.X # scaled data
    if scipy.sparse.issparse(X_exp):
        X_exp = X_exp.todense()
    LR = Ridge(fit_intercept=False,alpha=1.0)
    import pdb; pdb.set_trace()
    LR.fit(dummy,X_exp)

    if len(batch_key)>1:
        batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
    else:
        batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]
    
    dm = np.array(dummy)[:,batch_index]
    X_explained = dm.dot(LR.coef_[:,batch_index].T)
    X_remain = X_exp - X_explained
    ndata = sc.AnnData(X_remain)
    ndata.obs = adata.obs
    ndata.var = adata.var
    return ndata, X_explained

def regress_iter(adata,batch_key,confounder_key,bbknn_key,scale=True, approx = True,n_pcs = 50):
    if scale == True:
        sc.pp.scale(adata, max_value = 10)
    ndata, X_explained = regress_batch_v2(adata,batch_key=batch_key,confounder_key=confounder_key)
    sc.pp.pca(ndata)
    bbknn(ndata, batch_key = bbknn_key,n_pcs=n_pcs, approx=approx)
    return ndata, X_explained

start = time.time()
h5ad_obj = sc.read_h5ad('../3.results/3.anno_assess/data/merged.obj.UMAP.Scanpy.h5ad')

bdata, X_explained = regress_iter(h5ad_obj, ['percent.mt', 'percent.ribo', 'S.Score', 'G2M.Score'], ['Anno.Level.Fig.1'], 'orig.ident',scale = True)
sc.tl.umap(bdata)
bdata.write_h5ad("../3.results/3.anno_assess/data/merged.obj.UMAP.Scanpy.h5ad")
with plt.rc_context({'figure.figsize': (20, 12)}): 
    sc.pl.umap(bdata, color = 'Anno.Level.Fig.1', color_map = 'OrRd', size = 4, legend_loc = 'on data', show = False)
    plt.savefig(os.path.join('../3.results/3.anno_assess/figs', "Anno.UMAP.pdf"))
print(time.time() - start)

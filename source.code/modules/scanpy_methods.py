import os, sys
import scipy
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from sklearn.linear_model import Ridge

def diffusionMap(h5ad_file, n_pcs, meta_col, root_ident, root_idx = 0, use_rep = 'X_pca', out_dir = './'):
	# Read in h5ad object
	h5ad_obj = sc.read_h5ad(h5ad_file)

	# Compute a neighborhood graph of observations
	sc.pp.neighbors(h5ad_obj, n_neighbors = 40, n_pcs = int(n_pcs), use_rep = use_rep)

	# Set random naive cell as root for psuedotime
	h5ad_obj.uns['iroot'] = np.flatnonzero(h5ad_obj.obs[meta_col]  == root_ident)[int(root_idx)]

	# Compute diffusion map and pseudotime
	sc.tl.diffmap(h5ad_obj, n_comps = 3)
	sc.tl.dpt(h5ad_obj, n_branchings = 0, n_dcs = 3, allow_kendall_tau_shift = False, min_group_size = 10)

	# Export diffusion map coordinates and psuedotime values
	dat = h5ad_obj.obsm
	diffmap, diffmap_fil = dat['X_diffmap'], os.path.join(out_dir, "diffmap.txt")
	np.savetxt(diffmap_fil, diffmap, '%5.10f', delimiter = "\t")

	dpt, dpt_fil = h5ad_obj.obs.dpt_pseudotime, os.path.join(out_dir, "dpt.txt")
	dpt.to_csv(dpt_fil, sep = "\t", float_format = '%5.10f')
	return diffmap_fil, dpt_fil

def pagaAnalysis(h5ad_file, n_pcs, meta_col, use_rep = 'X_pca', res = 1.0):
	# Read in h5ad object
	h5ad_obj = sc.read_h5ad(h5ad_file)

	# Compute a neighborhood graph of observations
	sc.pp.neighbors(h5ad_obj, n_neighbors = 40, n_pcs = int(n_pcs), use_rep = use_rep)	
	
	# Clustering using louvain algorithm
	sc.tl.louvain(h5ad_obj, resolution = 1.0)

	# PAGA analysis
	sc.tl.paga(h5ad_obj, groups = 'louvain')
	sc.pl.paga(h5ad_obj, color = meta_col)
	return 0

def regress_batch_v2(exprs, meta_data, batch_key, confounder_key):
    '''batch regression tool
    batch_key=list of observation categories to be regressed out
    confounder_key=list of observation categories to be kept
    returns ndata with corrected X'''
    meta_data, confounder_key = pd.DataFrame(meta_data), [confounder_key]
    dummy = pd.get_dummies(meta_data[batch_key+confounder_key],drop_first=False)
    X_exp = exprs
    if scipy.sparse.issparse(X_exp):
        X_exp = X_exp.todense()
    LR = Ridge(fit_intercept=False,alpha=1.0)
    LR.fit(dummy,X_exp)

    if len(batch_key)>1:
        batch_index = np.logical_or.reduce(np.vstack([dummy.columns.str.startswith(x) for x in batch_key]))
    else:
        batch_index = np.vstack([dummy.columns.str.startswith(x) for x in batch_key])[0]

    dm = np.array(dummy)[:, batch_index]
    X_explained = dm.dot(LR.coef_[:, batch_index].T)
    X_remain = X_exp - X_explained
    return X_remain

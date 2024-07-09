import numpy as np
import scanpy as sc			   

#tdata = sc.read_h5ad('.obj.tcell.sub.UMAP.h5ad')
tdata = sc.read_h5ad('tcell.obj.UMAP.h5ad')
#tdata = sc.read_h5ad('.obj.t.diffumap.h5ad')
#tdata = sc.read_h5ad('./.obj.recluster.T.sub.h5ad') 
#import pdb; pdb.set_trace()
sc.tl.pca(tdata, svd_solver = 'arpack')
sc.pp.neighbors(tdata, n_pcs = 50, n_neighbors = 30)
#sc.pp.neighbors(tdata, n_pcs = 20, n_neighbors = 20)
sc.tl.umap(tdata)

sc.tl.diffmap(tdata)
tdata.uns['iroot'] = np.argmax(tdata.obsm['X_diffmap'][:, 1])
sc.tl.dpt(tdata)
tdata.obs.to_csv('.tcell.pseudo.xls', sep = '\t', header = True, index = True)

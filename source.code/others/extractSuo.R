import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_h5ad('Visium10X_data_TH.h5ad')
xx = adata[adata.obs['sample'] == 'WSSS_F_IMMsp9838711']
del xx.obsm['NMF']
del xx.obsm['means_cell_abundance_w_sf']
del xx.obsm['q05_cell_abundance_w_sf']
del xx.obsm['q95_cell_abundance_w_sf']
del xx.obsm['stds_cell_abundance_w_sf']

for key in list(xx.uns['spatial'].keys()):
	if key != 'spaceranger130_count_36654_WSSS_F_IMMsp9838711_GRCh38-2020-A':
		del xx.uns['spatial'][key]
		
xx.var_names = xx.var['SYMBOL'].astype(str)		
xx.var = xx.var[['feature_types', 'genome', 'mt']]

xx.write_h5ad("test.h5ad")
#Convert('test.h5ad', "Visium10X_data_TH.h5seurat", overwrite = TRUE)
#obj = LoadH5Seurat('Visium10X_data_TH.h5seurat', assays = 'data')
#count.data <- GetAssayData(obj) %>% as.data.frame
#colnames(count.data) <- gsub('spaceranger130_count_36654_WSSS_F_IMMsp9838711_GRCh38-2020-A_', '', colnames(count.data))
#obj <-  CreateSeuratObject(counts = count.data) %>% NormalizeData
#saveRDS(obj, file = 'obj_11.rds')

markers = ['EBI3', 'CCL17', 'CCR7', 'CSF2RB', 'CCL21', 'CCL22', 'TNFRSF18', 'CCL27', 'CXCL10', 'CXCL9', 'MS4A1', 'LAMP3']
mTEC = ['CCL19', 'KRT15', 'KRT19', 'FN1', 'SAA1', 'KRT14', 'PTGDS', 'CLU', 'KRT5', 'CH25H', 'KRT17', 'IGFBP5', 'CXCL14', 'TSC22D1', 'BST2', 'CALML3', 'IFITM3', 'CPE', 'APOE', 'AQP3'
]
anno = pd.read_csv('obj_11.xls', sep = '\t', header = 0)
anno.index = xx.obs.index
xx.obs['pred'] = anno['Pred']
sc.pl.spatial(xx, img_key="hires", color=["pred"], save = 'pred_obj_11.pdf', palette = ['#00B6EB', '#FB61D7'])
sc.pl.spatial(xx, img_key="hires", save = 'raw_obj_11.pdf')

import numpy as np
df = sc.get.obs_df(xx, 'CCL19')
ax = df.values
e = (ax - np.mean(ax)) / np.std(ax)
df = pd.DataFrame(df)
df.iloc[:, 0] = e
xx.obs['mTEC'] = df
sc.pl.spatial(xx, img_key="hires", color=["mTEC"], save = 'mTEC_obj_11.pdf', vmin = -1, vmax = 1)

#-------------------------------------

adata = sc.read_h5ad('Visium10X_data_TH.h5ad')
xx = adata[adata.obs['sample'] == 'WSSS_F_IMMsp9838716']
del xx.obsm['NMF']
del xx.obsm['means_cell_abundance_w_sf']
del xx.obsm['q05_cell_abundance_w_sf']
del xx.obsm['q95_cell_abundance_w_sf']
del xx.obsm['stds_cell_abundance_w_sf']

for key in list(xx.uns['spatial'].keys()):
	if key != 'spaceranger130_count_36811_WSSS_F_IMMsp9838716_GRCh38-2020-A':
		del xx.uns['spatial'][key]
xx.var_names = xx.var['SYMBOL'].astype(str)		
xx.var = xx.var[['feature_types', 'genome', 'mt']]

anno = pd.read_csv('obj_16.xls', sep = '\t', header = 0)
anno.index = xx.obs.index
xx.obs['pred'] = anno['Pred']
sc.pl.spatial(xx, img_key="hires", color=["pred"], save = 'pred_obj_16.pdf', palette = ['#00B6EB', '#FB61D7'])
sc.pl.spatial(xx, img_key="hires", save = 'raw_obj_16.pdf')

import numpy as np
df = sc.get.obs_df(xx, 'CCL19')
ax = df.values
e = (ax - np.mean(ax)) / np.std(ax)
df = pd.DataFrame(df)
df.iloc[:, 0] = e
xx.obs['mTEC'] = df
sc.pl.spatial(xx, img_key="hires", color=["mTEC"], save = 'mTEC_obj_16.pdf', vmin = -1, vmax = 1)


xx.write_h5ad("test.h5ad")
Convert('test.h5ad', "Visium10X_data_TH.h5seurat", overwrite = TRUE)
obj = LoadH5Seurat('Visium10X_data_TH.h5seurat', assays = 'data')
count.data <- GetAssayData(obj) %>% as.data.frame
colnames(count.data) <- gsub('spaceranger130_count_36811_WSSS_F_IMMsp9838716_GRCh38-2020-A_', '', colnames(count.data))
obj <-  CreateSeuratObject(counts = count.data) %>% NormalizeData
saveRDS(obj, file = 'obj_16.rds')
#---------------------------------------------------
adata = sc.read_h5ad('Visium10X_data_TH.h5ad')
xx = adata[adata.obs['sample'] == 'WSSS_F_IMMsp10864183']
del xx.obsm['NMF']
del xx.obsm['means_cell_abundance_w_sf']
del xx.obsm['q05_cell_abundance_w_sf']
del xx.obsm['q95_cell_abundance_w_sf']
del xx.obsm['stds_cell_abundance_w_sf']

for key in list(xx.uns['spatial'].keys()):
	if key != 'spaceranger130_count_WSSS_F_IMMsp10864183_GRCh38-2020-A':
		del xx.uns['spatial'][key]
xx.var_names = xx.var['SYMBOL'].astype(str)		
xx.var = xx.var[['feature_types', 'genome', 'mt']]

anno = pd.read_csv('obj_83.xls', sep = '\t', header = 0)
anno.index = xx.obs.index
xx.obs['pred'] = anno['Pred']
sc.pl.spatial(xx, img_key="hires", color=["pred"], save = 'pred_obj_83.pdf', palette = ['#00B6EB', '#FB61D7'])
sc.pl.spatial(xx, img_key="hires", save = 'raw_obj_83.pdf')

import numpy as np
df = sc.get.obs_df(xx, 'CCL19')
ax = df.values
e = (ax - np.mean(ax)) / np.std(ax)
df = pd.DataFrame(df)
df.iloc[:, 0] = e
xx.obs['mTEC'] = df
sc.pl.spatial(xx, img_key="hires", color=["mTEC"], save = 'mTEC_obj_83.pdf', vmin = -1, vmax = 1)

xx.write_h5ad("test.h5ad")

#Convert('test.h5ad', "Visium10X_data_TH.h5seurat", overwrite = TRUE)
#obj = LoadH5Seurat('Visium10X_data_TH.h5seurat', assays = 'data')
#count.data <- GetAssayData(obj) %>% as.data.frame
#colnames(count.data) <- gsub('spaceranger130_count_WSSS_F_IMMsp10864183_GRCh38-2020-A_', '', colnames(count.data))
#obj <-  CreateSeuratObject(counts = count.data) %>% NormalizeData
#saveRDS(obj, file = 'obj_83.rds')
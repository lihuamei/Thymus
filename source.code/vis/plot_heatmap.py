from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
import pdb

def draw_pseudo_heatmap_v2(adata,genelist,anno=None, clist=None, anno_list =None, gridspec_kw = {'height_ratios': [10,1]},figsize=(3,4),fontsize=12, mean_window = 100, anno_window=10):    
	import sklearn
	def running_mean(x, N):
		cumsum = np.cumsum(np.insert(x, 0, 0)) 
		return (cumsum[N:] - cumsum[:-N]) / float(N)
    
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.colors import LinearSegmentedColormap
	from matplotlib import colors
	fig,(axes) = plt.subplots(2, 1,  sharex= True,
                        gridspec_kw=gridspec_kw,
                        figsize=figsize)
	Xt = adata.X
	Rt = np.expm1(Xt)
	t = adata.obs['dpt_order']
	Rts = Rt[np.argsort(np.array(t))]
	bins = np.linspace(0,len(t),30)
	binlabel = np.digitize(np.array(t), bins)
	adata.obs['dpt_bins'] = [str(x) for x in binlabel]

	gidx = adata.var_names.isin(genelist)
	gene_order = list(adata.var_names[gidx])
	re_order = [gene_order.index(x) for x in genelist]
	
	normed = sklearn.preprocessing.normalize(Rts[:,gidx][:,re_order].todense().T,norm='max')
	runned = np.apply_along_axis(lambda x: running_mean(x,mean_window),1,normed)
	normed = sklearn.preprocessing.normalize(runned,norm='max')
	#import pdb; pdb.set_trace() 

	ax = axes[0]

	height=40
	im = ax.imshow(normed, extent = [0,20,0,height], vmax=1,origin='upper',cmap='viridis')
	#im = ax.imshow(normed, extent = [0,10,0,height], vmax=1,origin='upper',cmap='RdYlBu_r')
	ax.set_yticks(np.linspace(0,height,len(genelist)+1)[:-1]+height/(len(genelist)*2))
	ax.set_yticklabels(genelist[::-1],fontsize=fontsize)
	ax.set_xticks([])
	ax.grid(False)
	ax.set_aspect('auto')
	#fig.colorbar(im, cax=ax, orientation='horizontal') 
	if anno:

		anl = sorted(list(set(adata.obs[anno])))
		Ob = np.vstack([np.array(adata.obs[anno]==an) for an in anl]).astype(int).T
		Obs = Ob[np.argsort(np.array(t))]

		order = anl
		if anno_list!=None:
			re_order = [anl.index(x) for x in anno_list]
		else:
			re_order = [anl.index(x) for x in anl]
			anno_list = anl

		from sklearn.preprocessing import normalize
		normed = sklearn.preprocessing.normalize(Obs[:,re_order].T,norm='max')
		#runned = np.apply_along_axis(lambda x: running_mean(x,anno_window),1,normed)
		#normed = sklearn.preprocessing.normalize(runned,norm='max')
		normed = normed>0
		ax = axes[1]
        
		height=1.4


		#if clist:
		#	clist = ['white']+clist
		#else:
		#	import pdb; pdb.set_trace()
		#	clist = ['white']+list(np.array(adata.uns[anno+"_colors"])[re_order])

		for i in range(len(anl)):
			print(i, height/len(anl)*i,height/len(anl)*(i+1), anno_list[i])
			#import pdb; pdb.set_trace()
			#cmap = colors.ListedColormap([clist[0],clist[i+1]])
			try:
				cmap = colors.ListedColormap(['white', cols[anl[i]]])
			except:
				cmap = colors.ListedColormap(['white', 'blue'])
			#if i == 6: 
				#import pdb; pdb.set_trace()
			ax.imshow(normed[i:i+1], extent = [0,20,height/(len(anl))*i,height/(len(anl))*(i+1)], \
					vmax = 1, origin='upper', cmap=cmap, interpolation = 'none')
			# quadric
			#ax.imshow(normed[i:i+1], extent = [0,10,height/(len(anl))*i,height/(len(anl))*(i+1)], cmap=cmap)
		ax.set_yticks(np.linspace(0,height,len(anno_list)+1)[:-1]+height/(len(anno_list)*2))
		ax.set_yticklabels(anno_list,fontsize=fontsize)
		ax.set_xticks([])
		ax.set_ylim(0,height)
		ax.grid(False)
		ax.set_aspect('auto')

	plt.subplots_adjust(wspace=0, hspace=0.01)
    
	return im

def jp_save_fig(version,figcount,fig_format='pdf',fig_folder='11_Figs'):
    
    plt.savefig('%s/%s%s.%s'%(fig_folder,version,figcount,fig_format),bbox_inches='tight',format=fig_format,dpi=300)
    print('%s/%s%s.pdf'%(fig_folder,version,figcount))

import pyreadr
import scanpy as sc	
import numpy as np
tdata = sc.read_h5ad('TCR.tcell.obj.diffmap.h5ad')
fig3_anno_key = 'Anno.Level.11'

cols = pyreadr.read_r('.color.T.rds').get(None).to_dict()['Color']
T = tdata[tdata.obs[fig3_anno_key].isin(["DN_early", "DN_blast", "DN_re1", "DN_re2", "DN_ISP", "DP_blast3", "DP_blast2", "DP_blast1", "DP_re1", "DP_re2", 'abT(entry)', "CD8T", 'CD4T'])]

gene_list = ['CD34', 'IGLL1', 'TRGC2','TRDC','PTCRA','TRBC2','TRAC','CD4','CD8A','CD8B',
	'TCF4', 'NOTCH1', 'HES4',
	'CCR9','CCR7','CD5','CD27',
	'PCNA','CDK1','MKI67','TYMS','CCND2', 'CCND3','RAG1','RAG2','HRK','BMF','TP53INP1','ST18','HIVEP3','RGPD3','SMPD3','AQP3','RORC','CD2', 'SATB1','TOX2']
T.obs['dpt_order'] = np.argsort(np.argsort(T.obs['dpt_pseudotime']))

anno = 'Anno.Level.11'

anno_order = ["DN_early", "DN_blast", "DN_re1", "DN_re2", "DN_ISP", "DP_blast3", "DP_blast2", "DP_blast1", "DP_re1", "DP_re2", 'abT(entry)', "CD8T", 'CD4T']
anno_order = [x for x in anno_order if x in set(T.obs[anno])]
im = draw_pseudo_heatmap_v2(T,gene_list,anno=anno,anno_list = anno_order,figsize=(3,9),gridspec_kw={'height_ratios':[10,3]},fontsize=10,anno_window=3)
jp_save_fig('1.0','Fig2_H',fig_folder='./')

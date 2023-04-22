

print("checkpoint1", flush=True)
#######################################################################

DATA_DIR = '/lustre/groups/ml01/workspace/louis.kuemmerle/projects/A1/data2/' 
DATA_VERSION = 'oct22'
RESULTS_DIR = '/lustre/groups/ml01/workspace/louis.kuemmerle/projects/A1/results/'
SHAM = True
sham_str = '_wSham' if SHAM else ''

#######################################################################

import sys
sys.path.insert(0, "../../")

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
import numpy as np
import pandas as pd
import scanpy as sc
from util import pcr, upsample_underrepresented
print("checkpoint2", flush=True)
import itertools
from tqdm.notebook import tqdm
import json
from pathlib import Path


import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors

plt.rcParams['figure.dpi'] = 150
SMALL_SIZE = 17
MEDIUM_SIZE = 19
BIGGER_SIZE = 21
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#######################################################################

def pc_regressions(data,celltype,ct_key='cell_types',conditions="all",regions="all",uniform=False,n_comps=50,min_cells=20):
    """Run all pcr regressions we're interested in for one celltype
    
    uniform: bool
        Set this to False, first thought that would be important for weighting the samples. Instead this introduces a wrong
        "sample info weighting".
    """
    
    if conditions == "all":
        conditions = data.obs['condition'].unique().tolist()
    if regions == "all":
        regions = data.obs['region'].unique().tolist()
    
    results = {
        'full_set': dict({'condition':np.nan,'region':np.nan}, **{r:np.nan for r in regions}),
        'comb_subsets': {f'{comb[0]}_{comb[1]}':np.nan for comb in itertools.combinations(conditions,2)},
        'cond_subsets': {cond:{} for cond in conditions}
    }
    results_pval = {
        'full_set': dict({'condition':np.nan,'region':np.nan}, **{r:np.nan for r in regions}),
        'comb_subsets': {f'{comb[0]}_{comb[1]}':np.nan for comb in itertools.combinations(conditions,2)},
        'cond_subsets': {cond:{} for cond in conditions}
    }    
    
    iterations = (2+len(regions))
    iterations+= len([comb for comb in itertools.combinations(conditions,2)]) if (len(conditions) > 2) else 0
    iterations+= len(conditions) * (1+len(regions))
    pbar = tqdm(total=iterations)
        
    adata = data[data.obs[ct_key] == celltype].copy()
    
    # initialize covariates
    for r in regions:
        adata.obs[r] = 'other'
        adata.obs.loc[adata.obs['region'] == r,r] = r
    
    ###############################
    # regressions on full dataset #
    ###############################
    # - covariate = conditions (all conditions at once)
    # - covariate = regions:
    #    - all regions at once
    #    - one region vs all others (n_regions times)
    if not uniform:
        a = adata.copy()
        if a.n_obs > min_cells:
            if ('X_pca' in a.obsm): del a.obsm['X_pca']
            if ('PCs' in a.varm): del a.varm['PCs']
            if ('pca' in a.uns): del a.uns['pca']            
            sc.tl.pca(a, n_comps=min(n_comps,a.n_obs-1), svd_solver='arpack', use_highly_variable=False)
        
    if (a.n_obs > min_cells) and (a.obs['condition'].value_counts() > min_cells).all():
        if uniform:
            a = upsample_underrepresented(adata.copy(),'condition') 
        var_expl_and_pval = pcr(a,n_comps=n_comps,covariate='condition',verbose=False,recompute_pca=uniform)
        results['full_set']['condition'], results_pval['full_set']['condition'] = var_expl_and_pval
    else:
        results['full_set']['condition'] = np.nan
        results_pval['full_set']['condition'] = np.nan
    pbar.update()
    
    if (a.n_obs > min_cells) and (a.obs['region'].value_counts() > min_cells).all():
        if uniform:
            a = upsample_underrepresented(adata.copy(),'region')        
        var_expl_and_pval = pcr(a,n_comps=n_comps,covariate='region',verbose=False,recompute_pca=uniform)
        results['full_set']['region'], results_pval['full_set']['region'] = var_expl_and_pval
    else:
        results['full_set']['region'] = np.nan
        results_pval['full_set']['region'] = np.nan
    pbar.update()
    
    for r in regions:
        if (a.n_obs > min_cells) and (a.obs[r].value_counts() > min_cells).all():
            if uniform:
                a = upsample_underrepresented(adata.copy(),r)            
            var_expl_and_pval = pcr(a,n_comps=n_comps,covariate=r,verbose=False,recompute_pca=uniform)
            results['full_set'][r], results_pval['full_set'][r] = var_expl_and_pval
        else:
            results['full_set'][r] = np.nan
            results_pval['full_set'][r] = np.nan
        pbar.update()
    
    ######################################################
    # regress on each subset of 2 condition combinations #
    ###################################################### if number of conditions > 2
    if len(conditions) > 2:
        for comb in itertools.combinations(conditions, 2):
            a = adata[adata.obs['condition'].isin(comb)].copy()
            if (a.n_obs > min_cells) and (a.obs['condition'].value_counts() > min_cells).all():
                if uniform:
                    a = upsample_underrepresented(a,'condition')                
                var_expl_and_pval = pcr(a, n_comps=n_comps, covariate='condition', verbose=False, recompute_pca=True)
                results['comb_subsets'][f'{comb[0]}_{comb[1]}'], results_pval['comb_subsets'][f'{comb[0]}_{comb[1]}'] = var_expl_and_pval
            else:
                results['comb_subsets'][f'{comb[0]}_{comb[1]}'] = np.nan
                results_pval['comb_subsets'][f'{comb[0]}_{comb[1]}'] = np.nan
            pbar.update()
   
    ######################################
    # for each condition regress regions #
    ######################################
    for cond in conditions:
        a_cond = adata[adata.obs['condition'] == cond].copy()
        
        if not uniform:
            a = a_cond.copy()
            if a.n_obs > min_cells:
                if ('X_pca' in a.obsm): del a.obsm['X_pca']
                if ('PCs' in a.varm): del a.varm['PCs']
                if ('pca' in a.uns): del a.uns['pca']
                sc.tl.pca(a, n_comps=min(n_comps,a.n_obs-1), svd_solver='arpack', use_highly_variable=False)
        
        if (a.n_obs > min_cells) and (a.obs['region'].value_counts() > min_cells).all():
            if uniform:
                a = upsample_underrepresented(a_cond.copy(),'region')
            var_expl_and_pval = pcr(a,n_comps=n_comps,covariate='region',verbose=False,recompute_pca=uniform)
            results['cond_subsets'][cond]['region'], results_pval['cond_subsets'][cond]['region'] = var_expl_and_pval
        else:
            results['cond_subsets'][cond]['region'] = np.nan
            results_pval['cond_subsets'][cond]['region'] = np.nan
        pbar.update()
        
        for r in regions:
            if (a.n_obs > min_cells) and (a.obs[r].value_counts() > min_cells).all():            
                if uniform:
                    a = upsample_underrepresented(a_cond.copy(),r)
                var_expl_and_pval = pcr(a,n_comps=n_comps,covariate=r,verbose=False,recompute_pca=uniform)
                results['cond_subsets'][cond][r], results_pval['cond_subsets'][cond][r] = var_expl_and_pval
            else:
                results['cond_subsets'][cond][r] = np.nan
                results_pval['cond_subsets'][cond][r] = np.nan
            pbar.update()
        
    pbar.close()
    
    return results, results_pval


#######################################################################


def script(ct,ct_key):
    """
    
    ct_key: "level1" or "level2"
    """
    
    print("checkpoint3", flush=True)
    
    # Load adata
    adata = sc.read(DATA_DIR+f'cellxgene_{DATA_VERSION}{sham_str}_umaps.h5ad')
    adata = adata[adata.obs["region"] != "Scapula"]
    
    print("checkpoint4", flush=True)
    #########
    #cts = [ct for ct in adata.obs[ct_key].unique()]
    #cts = [ct for ct in adata.obs[ct_key].unique() if not os.path.isfile(results_dir+f'{ct}.json')]
    
    # THIS IS ACTUALLY FINE CODE; BUT NOT NEEDED RIGHT NOW!
    ## all regions (pool bones)
    #results_dir = RESULTS_DIR+f'pc_regression/{ct_key}/VS_bones_no_scapula/'
    #Path(results_dir).mkdir(parents=True, exist_ok=True)
    #cts = [ct for ct in adata.obs[ct_key].unique() if not os.path.isfile(results_dir+f'{ct}.json')]
    #a = adata.copy()
    #a.obs['region'] = a.obs['region'].astype(str)
    #a.obs.loc[~a.obs['region'].isin(['Brain','Meninges']),'region'] = 'Bones'    
    #for ct in tqdm(cts):
    #    if os.path.isfile(results_dir+f'{ct}.json') and os.path.isfile(results_dir+f'{ct}_pval.json'):
    #        continue
    #    results, results_pval = pc_regressions(a,ct,ct_key=ct_key,conditions="all",regions="all",uniform=False,n_comps=50,min_cells=20)
    #    with open(results_dir+f'{ct}.json', 'w') as file:
    #        json.dump(results, file)
    #    with open(results_dir+f'{ct}_pval.json', 'w') as file:
    #        json.dump(results_pval, file)
            
    # bones only
    results_dir = RESULTS_DIR+f'pc_regression/{ct_key}/bones_no_scapula/'
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    
    a = adata[~adata.obs['region'].isin(['Brain','Meninges'])].copy()
    
    if os.path.isfile(results_dir+f'{ct}.json') and os.path.isfile(results_dir+f'{ct}_pval.json'): 
        print("Results already exist.")
        return
    
    print("checkpoint5", flush=True)
    
    results, results_pval = pc_regressions(a,ct,ct_key=ct_key,conditions="all",regions="all",uniform=False,n_comps=50,min_cells=20)
    with open(results_dir+f'{ct}.json', 'w') as file:
        json.dump(results, file)
    with open(results_dir+f'{ct}_pval.json', 'w') as file:
        json.dump(results_pval, file)

    print("checkpoint6", flush=True)



#######################################################################

if __name__ == "__main__":
    
    print(sys.argv)
    
    ct_key = sys.argv[1]
    ct = sys.argv[2]
    ct = ct.replace("_"," ")
    
    script(ct,ct_key)

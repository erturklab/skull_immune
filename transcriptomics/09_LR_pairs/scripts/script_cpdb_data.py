import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from geosketch import gs

# Lists of celltypes
ct_selections = {
    'level1':['B cell', 'dendritic cell', 'erythroid cell', 'erythroid precursor', 'monocyte', 
              'neutrophil', 'NK cell', 'NK-T cell', 'progenitors', 'T cell'],
    'level2_min10':['Cd8 T cell', 'erythroblast', 'erythrocyte', 'erythroid cell', 
                    'erythroid progenitor', 'immature neutrophil', 'mature B cell', 
                    'mature neutrophil', 'monocyte progenitor', 'monocyte-primed GMP', 
                    'neutrophil-primed GMP', 'NK cell', 'NK-T cell', 'plasmacytoid DC', 
                    'pre neutrophil', 'pro neutrophil'],
    'level2_min5':['basophil', 'Cd4 T cell', 'Cd8 T cell', 'common myeloid progenitor', 'erythroblast', 
                   'erythrocyte', 'erythroid cell', 'erythroid progenitor', 'granulocyte-monocyte progenitor', 
                   'immature B cell', 'immature neutrophil', 'mature B cell', 'mature neutrophil', 
                   'monocyte progenitor', 'monocyte-primed GMP', 'neutrophil-primed GMP', 'NK cell', 'NK-T cell', 
                   'non-classical monocyte', 'plasmacytoid DC', 'pre neutrophil', 'pro neutrophil'],
}

def sample_n_per_celltype(adata, n=1000, ct_key="level1", seed=0):
    """
    
    gemoetric sketching is only used for seed == 0, otherwise random sampling.
    
    """
    obs_names = []
    
    n = min(n, adata.obs[ct_key].value_counts().max())
    
    for ct in adata.obs[ct_key].unique():
        ct_filt = (adata.obs[ct_key] == ct)
        adata_ct = adata[ct_filt, :].copy()
        if adata_ct.n_obs < n:
            # upsampling
            # douplicate (repeat_n times) observerations if n_obs*2 < n
            repeat_n = n // adata_ct.n_obs
            obs_names += adata_ct.obs_names.tolist() * repeat_n
            # sample remaining number of observations
            sample_n_obs = n % adata_ct.n_obs
            np.random.seed(seed)
            adata_ct = adata_ct[np.random.choice(adata_ct.n_obs, sample_n_obs, replace=False), :]
        elif (seed==0):
            # downsampling with gemoetric sketching
            sc.tl.pca(adata_ct)
            X_dimred = adata_ct.obsm['X_pca']
            sketch_index = gs(X_dimred, n, replace=False,)
            adata_ct = adata_ct[sketch_index, :]
        else:
            # downsampling with random sampling
            np.random.seed(seed)
            adata_ct = adata_ct[np.random.choice(adata_ct.n_obs, n, replace=False), :]
        obs_names += adata_ct.obs_names.tolist()

    return obs_names
            


    
def create_data_txt(condition="all",
                    region="all",
                    ct_selection_key='level1',
                    test=0,
                    data_dir='./data',
                    dataset="cellxgene_april21_wSham_umaps",
                    seed=0,
                    subsample_per_ct=400,
                    ):
    """
    
    """
    # Define subset of celltypes used for interactions
    CELLTYPES = ct_selections[ct_selection_key]

    if ct_selection_key[:6] == 'level1':
        CT_KEY = "level1"
    elif ct_selection_key[:6] == 'level2':
        CT_KEY = "level2"
    else:
        CT_KEY = "cell_types"

    print("### Check: Included celltypes ###")
    print(f"CT_KEY: {CT_KEY}")
    print("CELLTYPES:")
    print(CELLTYPES)

    if test == 1:
        dest_bp = f"{data_dir}/cellphonedb/{ct_selection_key}_{condition}_{region}_seed{seed}_test/"
    elif test == 0:
        dest_bp = f"{data_dir}/cellphonedb/{ct_selection_key}_{condition}_{region}_seed{seed}/"
    dest_data_csv = dest_bp + f"{dataset}_non_log.txt"
    dest_meta_csv = dest_bp + f"{dataset}_non_log_meta.txt"

    if (not (os.path.isfile(dest_data_csv))) and (not (os.path.isfile(dest_meta_csv))):

        data_path = f"{data_dir}/{dataset}.h5ad"
        adata = sc.read(data_path)
        adata = adata[np.in1d(adata.obs[CT_KEY], CELLTYPES), :]
        print(adata)
        
        if test == 1:
            a = adata[::20, :].copy()
        elif test == 0:
            a = adata.copy()
        a.var.index = a.var.index.astype(str).str.upper()
        print(a)
        if condition != 'all':
            cond_filt = (a.obs['condition'] == condition)
            a = a[cond_filt, :]
        print(a)
        if region != 'all':
            region_filt = (a.obs['region'] == region)
            a = a[region_filt, :]
            
        obs_sample = sample_n_per_celltype(a, n=subsample_per_ct, ct_key=CT_KEY, seed=seed)
        a = a[obs_sample,:]
            
        print(a)
        print("### Subsetting Check ###")
        print(f"Specified condition: {condition}")
        print(f"\tUnique conditions: {a.obs['condition'].unique()}")
        print(f"Specified region: {region}")
        print(f"\tUnique regions: {a.obs['region'].unique()}")

        #a.X = np.expm1(a.X) 
        df_expr_matrix = a.X.copy()
        df_expr_matrix = df_expr_matrix.T
        df_expr_matrix = pd.DataFrame(df_expr_matrix.toarray())
        df_expr_matrix.columns = a.obs.index
        df_expr_matrix.set_index(a.var.index, inplace=True)
        Path(dest_bp).mkdir(parents=True, exist_ok=True)
        df_expr_matrix.to_csv(dest_data_csv, sep='\t')
        df_meta = pd.DataFrame(
            data={'Cell': list(a.obs.index), 'cell_type': list(a.obs[CT_KEY])})
        df_meta.set_index('Cell', inplace=True)
        df_meta.to_csv(dest_meta_csv, sep='\t')
    else:
        print("Data files")
        print(f"\t{dest_data_csv}")
        print(f"\t{dest_meta_csv}")
        print("already exist.")


if __name__ == "__main__":
    """
    
    """
    data_dir = str(sys.argv[1])
    dataset = str(sys.argv[2])

    ct_selection_key = str(sys.argv[3])
    condition = str(sys.argv[4])
    region = str(sys.argv[5])

    # test run (either 0 - no, or 1 - yes)
    test = int(sys.argv[6])
    
    seed = int(sys.argv[7])

    print(f"{data_dir=}, {dataset=}, {ct_selection_key=}, {condition=}, {region=}, {test=}, {seed=}")

    create_data_txt(
        condition=condition,
        region=region,
        ct_selection_key=ct_selection_key,
        test=test,
        data_dir=data_dir,
        dataset=dataset,
        seed=seed,
    )

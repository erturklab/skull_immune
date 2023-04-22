import numpy as np
import pandas as pd
import scanpy as sc
import sklearn
from scipy import sparse
from scipy.sparse import issparse



##############
# Data utils #
##############

def upsample_underrepresented(adata,covariate):
    """Upsample adata observations to have equal occurances of the given covariate values
    
    Thoroughly tested.
    
    adata: AnnData
    covariate: str
        Categorical adata.obs key.
    """
    counts_df = adata.obs[covariate].value_counts()
    ref_val = counts_df.idxmax()
    ref_count = counts_df.loc[ref_val]
    vals = [v for v in counts_df.index if v != ref_val]
    counts = counts_df.loc[vals].values.tolist()
    
    for val,count in zip(vals,counts):
        n_duplicates = (ref_count // count) - 1
        n_rest = ref_count % count
        obs_filter = (adata.obs[covariate] == val)
        adata = adata.concatenate([adata[obs_filter] for _ in range(n_duplicates)] + [adata[obs_filter][:n_rest]])
    return adata



#################
# PC Regression #
################# Code from scIB.metrics (we cut out stuff that we don't need)

def pcr(adata, covariate, n_comps=50, recompute_pca=True, verbose=False, n_perm=1000):
    """
    PCR for Adata object
    Checks whether to
        + compute PCA on expression data
        + use existing PCA (only if PCA entry exists)
        + recompute PCA on expression matrix (default)
    params:
        adata: Anndata object
        n_comps: number of PCs if PCA should be computed
        covariate: key for adata.obs column to regress against
    return:
        R2Var of PCR
    """
        
    if verbose:
        print(f"covariate: {covariate}")
    batch = adata.obs[covariate]
        
    # use existing PCA computation
    if (recompute_pca == False) and ('X_pca' in adata.obsm) and ('pca' in adata.uns):
        if verbose:
            print("using existing PCA")
        return pc_regression(adata.obsm['X_pca'], batch, pca_var=adata.uns['pca']['variance'], n_perm=n_perm)
    
    # recompute PCA
    else:
        if verbose:
            print(f"compute PCA n_comps: {n_comps}")
        return pc_regression(adata.X, batch, n_comps=n_comps, n_perm=n_perm)

def pc_regression(data, variable, pca_var=None, n_comps=50, svd_solver='arpack', verbose=False, n_perm=1000, seed=0):
    """
    params:
        data: expression or PCA matrix. Will be assumed to be PCA values, if pca_sd is given
        variable: series or list of batch assignments
        n_comps: number of PCA components for computing PCA, only when pca_sd is not given. If no pca_sd is given and n_comps=None, comute PCA and don't reduce data
        pca_var: iterable of variances for `n_comps` components. If `pca_sd` is not `None`, it is assumed that the matrix contains PCA values, else PCA is computed
    PCA is only computed, if variance contribution is not given (pca_sd).
    """

    if isinstance(data, (np.ndarray, sparse.csr_matrix)):
        matrix = data
    else:
        raise TypeError(f'invalid type: {data.__class__} is not a numpy array or sparse matrix')

    # perform PCA if no variance contributions are given
    if pca_var is None:

        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = 'full'
            if issparse(matrix):
                matrix = matrix.toarray()

        if verbose:
            print("compute PCA")
        pca = sc.tl.pca(matrix, n_comps=n_comps, use_highly_variable=False,
                        return_info=True, svd_solver=svd_solver, copy=True)
        X_pca = pca[0].copy()
        pca_var = pca[3].copy()
        del pca
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]

    ## PC Regression
    if verbose:
        print("fit regression on PCs")

    # handle categorical values
    if pd.api.types.is_numeric_dtype(variable):
        variable = np.array(variable).reshape(-1, 1)
    else:
        if verbose:
            print("one-hot encode categorical values")
        variable = pd.get_dummies(variable)

    # fit linear model for n_comps PCs
    r2 = []
    for i in range(n_comps):
        pc = X_pca[:, [i]]
        lm = sklearn.linear_model.LinearRegression()
        lm.fit(variable, pc)
        r2_score = np.maximum(0,lm.score(variable, pc))
        r2.append(r2_score)

    Var = pca_var / sum(pca_var) * 100
    R2Var = sum(r2*Var)/100
    
    
    ### permutation test
    rng = np.random.default_rng(seed)
    pval = None
    if n_perm:
        R2Var_perm = []
        
        # measure variance explained over n permutations
        for p in range(n_perm):

            # permute covariate
            perm = rng.permutation(len(variable))
            if isinstance(variable,pd.DataFrame):
                variable = variable.iloc[perm]
            else:
                variable = variable[perm]
                
            # regression
            r2 = []
            for i in range(n_comps):
                pc = X_pca[:, [i]]
                lm = sklearn.linear_model.LinearRegression()
                lm.fit(variable, pc)
                r2_score = np.maximum(0,lm.score(variable, pc))
                r2.append(r2_score)

            Var = pca_var / sum(pca_var) * 100
            R2Var_perm.append(sum(r2*Var)/100)
        
        pval = np.sum((np.array(R2Var_perm) > R2Var))/n_perm
        
        #import seaborn as sns
        #import matplotlib.pyplot as plt
        #
        #sns.distplot(R2Var_perm,kde=False)
        #plt.axvline(R2Var, color='k', linestyle='dashed', linewidth=1)
        #plt.show()

    return R2Var, pval



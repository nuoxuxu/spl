import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scquint.data import calculate_PSI
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns

# First case: clustering cells at the t-type level, keeping only introns that are found in all t-types

adata = anndata.read_h5ad("/nethome/kcni/nxu/spl/results/markers_spl.h5ad")
ephys_data = pd.read_csv("/nethome/kcni/nxu/KCNC1/pseudobulk_by_ttypes.csv", header = 0, index_col = 0)
ephys_feats_list = ['avg_rate', 'width','threshold_v', 'peak_v', 'input_resistance', 'fi_fit_slope', 'sag', 'n','v_baseline']
ion_channel_genes = pd.read_html("https://www.guidetopharmacology.org/GRAC/IonChannelListForward?class=VGIC", header = 1)[0]['Mouse gene symbol'].unique()

# keeping cell types that have corrersponding ephys values
adata = adata[adata.obs.index.isin(ephys_data.markers.values)]
# keeping introns that are associated with ion channel genes
adata = adata[:, adata.var.gene_name.isin(ion_channel_genes)] # type: ignore
adata.layers["PSI_raw"] = calculate_PSI(adata)
#calculate the global mean and standard deviations for PSI_raw
global_mean = np.nanmean((np.nanmean(adata.layers["PSI_raw"], axis = 0)))
global_std = np.nanmean((np.nanstd(adata.layers["PSI_raw"], axis = 0))  # type: ignore
# filtering for introns that have higher mean and std than global values
adata = adata[:, np.logical_and(np.greater(np.nanstd(adata.layers["PSI_raw"], axis = 0), global_std), np.greater(np.nanmean(adata.layers["PSI_raw"], axis = 0), global_mean))]
# keeping only introns that are expressed in all t-types
ttypes_by_introns = pd.DataFrame(adata.layers["PSI_raw"], index = adata.obs_names, columns = adata.var_names).dropna(axis = 1)

intron_by_ephys_coor = pd.DataFrame(
        {ephys_feat: np.apply_along_axis(lambda x: pearsonr(x, ephys_data[ephys_feat]), axis = 0, arr = ttypes_by_introns)[0] for ephys_feat in ephys_feats_list}, 
        index = adata[:, ~np.isnan(adata.layers["PSI_raw"]).any(axis=0)].var_names)

intron_by_ephys_pvalue = pd.DataFrame(
        {ephys_feat: np.apply_along_axis(lambda x: pearsonr(x, ephys_data[ephys_feat]), axis = 0, arr = ttypes_by_introns)[1] for ephys_feat in ephys_feats_list}, 
        index = adata[:, ~np.isnan(adata.layers["PSI_raw"]).any(axis=0)].var_names)

sns.heatmap(
        intron_by_ephys_coor, 
        annot = np.vectorize({True: "*", False: " "}.get)(fdrcorrection(intron_by_ephys_pvalue.to_numpy().flatten(order = "C"))[0].reshape(intron_by_ephys_pvalue.shape)), 
        fmt = '')

# features associated with selected introns

data = adata.var.assign(sig_coor = adata.var.index.isin(ttypes_by_introns.columns))[["total_unique_mapping", "color"]]
data["total_unique_mapping"] = np.log(data["total_unique_mapping"])
sns.stripplot(data = data, x = "total_unique_mapping", hue = "color")

#-------------------------------------------------------------
# Second case: raw intron counts at the cell-level, smooth PSI

adata = anndata.read_h5ad("/external/rprshnas01/netdata_kcni/stlab/Nuo/output/SJ/raw/adata_spl.h5ad")
adata = adata[:, adata.var.gene_name.isin(ion_channel_genes)] # type: ignore
adata.layers["PSI_smooth"] = calculate_PSI(adata, smooth = True)

cells_by_introns = pd.DataFrame(adata.layers["PSI_smooth"], index = adata.obs_names, columns = adata.var_names)

ephys_data = pd.read_csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/ephys/processed/mouse_ephys_features_gouwens.csv")

intron_by_ephys_coor = pd.DataFrame(
        {ephys_feat: np.apply_along_axis(lambda x: pearsonr(x, ephys_data[ephys_feat]), axis = 0, arr = cells_by_introns)[0] for ephys_feat in ephys_feats_list}, 
        index = adata[:, ~np.isnan(adata.layers["PSI_smooth"]).any(axis=0)].var_names)

intron_by_ephys_pvalue = pd.DataFrame(
        {ephys_feat: np.apply_along_axis(lambda x: pearsonr(x, ephys_data[ephys_feat]), axis = 0, arr = cells_by_introns)[1] for ephys_feat in ephys_feats_list}, 
        index = adata[:, ~np.isnan(adata.layers["PSI_smooth"]).any(axis=0)].var_names)


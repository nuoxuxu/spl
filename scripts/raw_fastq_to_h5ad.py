from pathlib import Path
from src.helper_functions import add_metadata_to_barcodes, create_h5ad

output_dir = Path("/external/rprshnas01/netdata_kcni/stlab/Nuo/output/SJ/raw")
with open("../results/corrupted_gz_removed.txt", "r") as f:
    lines = f.readlines()
directory_list = [line.strip('\n') for line in lines]

# transform the directory list to a filemanifest.tsv
file_manifest = "\n".join(
    [(lambda directory: "\t".join([str(file) for file in Path(directory).iterdir()][::-1] +
    [Path(directory).name]))(directory_path) for directory_path in directory_list]
    )
with open("../results/file_manifest_STARsolo.tsv", 'w') as f:
    f.write(file_manifest)

#SJ
add_metadata_to_barcodes(directory_list, output_dir)

create_h5ad(output_dir, "barcodes_ann.tsv")

#Gene
add_metadata_to_barcodes(directory_list, "/external/rprshnas01/netdata_kcni/stlab/Nuo/output/Gene/raw")

#-----------------------------------------------------------------------

#TODO remove cell barcodes that do not have corresponding ephys features
from scquint.data import load_adata_from_starsolo, add_gene_annotation, group_introns
import pandas as pd
import numpy as np
adata = load_adata_from_starsolo(
    Path("/external/rprshnas01/netdata_kcni/stlab/Nuo/output/SJ/raw"), obs_filename="barcodes_ann.tsv")
ephys_data = pd.read_csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/ephys/processed/mouse_ephys_features_gouwens.csv")
file_manifest = pd.read_csv(
    "/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/2021-09-13_mouse_file_manifest.csv",
    dtype = {"cell_specimen_id" : str}
    )
temp = pd.merge(
    pd.merge(ephys_data, file_manifest[["cell_specimen_id", "file_name"]], how = "left", left_on = "ephys_session_id", right_on = "file_name"),
    file_manifest.loc[file_manifest.file_type == "reverse_fastq", ["cell_specimen_id", "file_name"]],
    how = "left", 
    on = "cell_specimen_id",
)
temp = temp.assign(barcode = temp.file_name_y.apply(lambda barcode : str.rsplit(barcode, "_", 1)[0]))
adata = adata[adata.obs.index.isin(temp.barcode.values)]

adata = add_gene_annotation(adata, "/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Ensembl/Mouse/Release_104/Raw/Mus_musculus.GRCm39.104.gtf")
adata = group_introns(adata, by="three_prime")
adata.write_h5ad("/nethome/kcni/nxu/spl/results/raw_spl_unlabelled_removed.h5ad")



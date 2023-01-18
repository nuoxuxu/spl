# Importing packages
from pathlib import Path
from scquint.data import load_adata_from_starsolo, add_gene_annotation, group_introns
import pandas as pd
import hashlib
import numpy as np
import scipy

def add_metadata_to_barcodes(directory_list, output_dir):
    """Write a tsv file that contains metadata info for cell samples

    Parameter
    ---------
    directory_list : list
        a list of strings that identify each cell samples, e.g., 'PS0810_E1-50_GCTCATGA-TCTCTCCG',
        either generated by get_directory_list.sh (prefix has to be removed),
        or read from barcodes.tsv in STARsolo output
    outputdir : str
        Path-like string to the output directory

    Output
    -------
    barcodes_ann.tsv
        a tsv file that contains the identifying strings as file_name column, and extra columns of metadata annotation
    """
    metadata = pd.read_csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/20200711_patchseq_metadata_mouse.csv")
    mouse_file_manifest = pd.read_csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/2021-09-13_mouse_file_manifest.csv")
    mouse_file_manifest = mouse_file_manifest.loc[(mouse_file_manifest['technique'] == 'transcriptomics') & ((mouse_file_manifest['file_type'] == 'reverse_fastq') | (mouse_file_manifest['file_type'] == 'forward_fastq'))]
    mouse_file_manifest['cell_specimen_id'] = mouse_file_manifest['cell_specimen_id'].astype('int64')
    mouse_file_manifest['file_name'] = mouse_file_manifest['file_name'].apply(lambda x : x.rsplit("_", 1)[0])
    barcodes = pd.DataFrame([Path(directory).name for directory in directory_list], columns = ['file_name'])
    new_barcodes = pd.merge(left = barcodes, right = mouse_file_manifest[['file_name', 'cell_specimen_id']], how = 'left', on = 'file_name')
    new_barcodes = new_barcodes.drop_duplicates(subset = ['file_name'], ignore_index = True)
    new_barcodes = pd.merge(left = new_barcodes, right = metadata[['cell_specimen_id', 'Unnamed: 21']], how = 'left', on = 'cell_specimen_id')
    new_barcodes = new_barcodes.rename(columns = {'Unnamed: 21': 'markers'})
    new_barcodes['markers'] = new_barcodes['markers'].astype(str)
    new_barcodes['subclass'] = new_barcodes['markers'].map(lambda x : x.split()[0])
    new_barcodes.to_csv(Path(output_dir).joinpath("barcodes_ann.tsv"), sep = "\t", index = False)

def write_features_and_mtx(directory_list, output_dir, suffix):
    """Write features.tsv and matrix.mtx from a list of paths to cell directories containing gene or SJ counts

    Parameter
    ---------
    directory_list : list
        A list of absolute paths to directories containing gene or SJ counts
    
    output_dir : str
        Path-like string to the output directory

    suffix : str
        Either '_SJ.out.tab' for SJs, or '_ReadsPerGene.out.tab' for gene counts
    
    Output
    ------
    Write features.tsv and matrix.mtx in the output directory
    """
    if suffix == "_SJ.out.tab":
        list_of_df = []
        for i, directory in zip(range(len(directory_list)), directory_list):
            tab_path = Path(directory).joinpath("".join([Path(directory).name, suffix]))
            tab = pd.read_csv(tab_path, sep="\t", names = ["chromosome", "start", "end", "strand", "intron_motif", "annotation", "unique", "multi", "max_overhang"])
            tab = tab.assign(barcodes = i)
            list_of_df.append(tab)
            print(f"{directory}added")
        combined = pd.concat(list_of_df, ignore_index=True)
        print("All SJ.out.tab files combined")
        features = combined.drop(columns = 'barcodes').copy()
        features = features.drop_duplicates(subset=["chromosome", "start", "end", "intron_motif", "strand", "annotation"])
        features = pd.merge(features, combined.groupby(['start', 'end']).aggregate({'unique':np.sum, 'multi':np.sum, 'max_overhang':max}), how = 'inner', on = ['start', 'end'])
        features = features.drop(columns = ['unique_x', 'multi_x', 'max_overhang_x'])
        features = features.rename(columns = {'unique_y':'unique', 'multi_y':'multi', 'max_overhang_y':'max_overhang'})

        with open(Path(output_dir).joinpath("features.tsv"), 'w') as f:
            features.to_csv(f, sep = '\t', header = None, index = False)
        print("features.tsv written")
        combined = pd.merge(combined, features[['start', 'end']].assign(row = features.index), how = 'inner', on = ['start', 'end'])
        mtx = scipy.sparse.coo_matrix((combined['unique'].to_numpy(), (combined['row'].to_numpy(), combined['barcodes'].to_numpy())), dtype=np.int64)
        with open(Path(output_dir).joinpath("matrix.mtx"), "wb") as f:
            scipy.io.mmwrite(f, mtx)
        print("matrix.mtx written")
    if suffix == "_ReadsPerGene.out.tab":
        list_of_df = []
        for i, directory in zip(range(len(directory_list)), directory_list):
            tab_path = Path(directory).joinpath("".join([Path(directory).name, suffix]))
            tab = pd.read_csv(tab_path, sep="\t", usecols=[1], skiprows=[0, 1, 2, 3], names=['unique'])
            tab = tab.assign(barcodes = i)
            tab = tab.assign(row = tab.index)
            list_of_df.append(tab)
            print(f"{directory}added")
        combined = pd.concat(list_of_df, ignore_index=True)
        print("All ReadsPerGene.out.tab files combined")
        mtx = scipy.sparse.coo_matrix((combined['unique'].to_numpy(), (combined['row'].to_numpy(), combined['barcodes'].to_numpy())), dtype=np.int64)
        with open(Path(output_dir).joinpath("matrix.mtx"), "wb") as f:
            scipy.io.mmwrite(f, mtx)
        print("matrix.mtx written")

def create_h5ad(target_dir, obs_filename):
    """Create h5ad file from the STARsolo output

    Parameter
    ---------
    target_dir :  str
        The path to the directory that contains the three files (features.tsv, barcodes.tsv, matrix.mtx)
    obs_filename : str
        Default: barcodes.tsv, could be barcodes_ann.tsv with the annotations

    Output
    ------
    Write a h5ad file: adata_spl.h5ad
    """
    adata = load_adata_from_starsolo(Path(target_dir), obs_filename=obs_filename)
    adata = add_gene_annotation(adata, "/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Ensembl/Mouse/Release_104/Raw/Mus_musculus.GRCm39.104.gtf")
    adata = group_introns(adata, by="three_prime")
    adata.write_h5ad(Path(target_dir).joinpath('adata_spl.h5ad'))

def remove_corrupted_gz(directory_list):
    """Remove cell directories that contain corrupted .gz files

    Parameter
    ---------
    directory_list : list
        A list of paths to the cell directories that contain raw fastq files

    Output
    ------
    new_directory_list : list
        A list of paths to the directories, without directories that contain corrupted .gz files
    """

    file_manifest = pd.read_csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/2021-09-13_mouse_file_manifest.csv")
    new_directory_list = []
    for directory in directory_list:
        truth_values = np.array([])
        for filename in Path(directory).iterdir():
            hash_from_manifest = file_manifest.loc[file_manifest['file_name'] == filename.name, 'md5_checksum']
            with open(filename, 'rb') as f:
                bytes = f.read()
                calculated_hash = hashlib.md5(bytes).hexdigest()
                truth_values = np.append(truth_values, hash_from_manifest == calculated_hash)
        if np.all(truth_values):
            new_directory_list.append(directory)
    return new_directory_list
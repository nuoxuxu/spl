STAR --runThreadN 12 \
--genomeDir /external/rprshnas01/netdata_kcni/stlab/Genomic_references/Ensembl/Mouse/Release_104/USE_THIS_genomeDir/ \
--readFilesCommand zcat \
--soloType SmartSeq \
--soloUMIdedup Exact \
--readFilesManifest /nethome/kcni/nxu/scQuint/new_manifest.tsv \
--soloOutFileNames all_cells/ features.tsv barcodes.tsv matrix.mtx \
--soloStrand Unstranded \
--soloFeatures Gene SJ




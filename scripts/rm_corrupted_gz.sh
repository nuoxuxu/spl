cd /external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/transcriptomics/raw
arr=()
for directory in $(ls -d */ | head -10); do
    cd $directory
    for file in *gz; do
        if gzip -t $file; then
            echo "this $file file is not corrupted"
        fi
        arr+=$file
    done
    cd ..
done
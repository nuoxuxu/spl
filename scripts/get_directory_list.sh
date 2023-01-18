#!/usr/bin/env bash
ls -d /external/rprshnas01/netdata_kcni/stlab/Nuo/patchseqBams/STAR_results/*/ \
| egrep '[ATGC]{8,}/$' \
| sed -E -e '/SM-GBMG3_E1-50_AGTTAGCTGG-CATTCTCATC\/$/d' - \
> $(pwd)/../results/directory_list.txt
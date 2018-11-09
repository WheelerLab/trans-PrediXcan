#!/bin/bash
#PBS -N cat-mtpx
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

###concatenate output from 65_multi_trans-px_allgenes-PCs.R 

rm o p q

zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_overall_results_list000_2017-11-21.txt |head -n 1 > allheader

for (( i = 0 ; i <= 9; i++))
do
    echo $i
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_overall_results_list00${i}_2017-*.txt.gz | tail -n +2 | grep -v y >>o
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_DGN_overall_results_list00${i}_2017-*.txt.gz | tail	-n +2 >>p
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_GEU_overall_results_list00${i}_2017-*.txt.gz | tail	-n +2 >>q
done

for (( i = 10 ; i <= 99; i++))
do
    echo $i
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_overall_results_list0${i}_2017-*.txt.gz | tail -n +2 | grep -v y >>o
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_DGN_overall_results_list0${i}_2017-*.txt.gz | tail  -n +2 >>p
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_GEU_overall_results_list0${i}_2017-*.txt.gz | tail  -n +2 >>q
done

for (( i = 100 ; i <= 300; i++))
do
    echo $i
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_overall_results_list${i}_2017-*.txt.gz | tail -n +2 | grep -v y >>o
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_DGN_overall_results_list${i}_2017-*.txt.gz | tail  -n +2 >>p
    zcat multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_GEU_overall_results_list${i}_2017-*.txt.gz | tail  -n +2 >>q
done
    
cat allheader o > multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_overall_results_2017-12-11.txt
cat allheader p	> multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_DGN_overall_results_2017-12-11.txt
cat allheader q	> multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_GEU_overall_results_2017-12-11.txt

gzip multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_*_overall_results_2017-12-11.txt

rm o p q allheader

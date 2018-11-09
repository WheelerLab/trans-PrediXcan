#!/usr/bin/env python

'''make a run script for each chr and output a qsub file'''

qsubfile = open('../qsub.txt','w')
prescript = '65_multi_trans-px_allgenes-PCs'
pop = 'GEU'

for i in range(301):
    num = str(i).zfill(3)
    outfilename = 'run_' + prescript + '_' + num + '.sh'
    outfile = open(outfilename,'w')
    output = '''#!/bin/bash
#PBS -N mt-tpx.''' + pop + num +'''\n#PBS -S /bin/bash
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load R/3.3.2

time R --no-save < ''' + prescript + '.R --args ' + num + ' ' + pop + '\n'
    outfile.write(output)
    qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')



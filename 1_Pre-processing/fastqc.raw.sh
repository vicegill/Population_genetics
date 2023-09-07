#!/bin/bash
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=4:00:00
#SBATCH -J fastqc.raw.data


for i in /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/*fq.gz; 
do fastqc $i;
done


multiqc .


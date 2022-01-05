#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --job-name=ca_r1
#SBATCH --mem=18G
#SBATCH --cpus-per-task=4
#SBATCH --output=/mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/outputs/ca_r1.out

echo "####Started /mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/scripts/07ca_r1.sh at $(date)" >> outputs/log.txt
cd /mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1 || exit

## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.

mkdir -p /mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/outputs/cutadapt && \
 less r1.fastq.xz | cutadapt -  \
    -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC  -a GTCATAGCTGTTT \
  -e 0.1 -n 3 \
  -m 8 -M 42 \
  --too-short-output=/mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/outputs/cutadapt/r1_tooshort.fastq \
  --too-long-output=/mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/outputs/cutadapt/r1_toolong.fastq \
  --untrimmed-output=/mnt/sshfs/cbcbsub01/fs/cbcb-lab/nelsayed/scratch/atb/dnaseq/paeruginosa_stamp_2020/preprocessing/a1/outputs/cutadapt/r1_untrimmed.fastq \
  2>outputs/cutadapt.err 1>r1-trimmed_ca.fastq

## The following lines give status codes and some logging
echo $? > outputs/status/ca_r1.status
echo "###Finished ${SLURM_JOBID} 07ca_r1.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

##walltime=$(qstat -f -t "${SLURM_JOBID}" | grep 'resources_used.walltime' | awk -F ' = ' '{print $2}')
##echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
##mem=$(qstat -f -t | grep "${SLURM_JOBID}" | grep 'resources_used.mem' | awk -F ' = ' '{print $2}')
##echo "#### memory used by ${SLURM_JOBID} was: ${mem:-null}" >> outputs/log.txt
##vmmemory=$(qstat -f -t "${SLURM_JOBID}" | grep 'resources_used.vmem' | awk -F ' = ' '{print $2}')
##echo "#### vmemory used by ${SLURM_JOBID} was: ${vmmemory:-null}" >> outputs/log.txt
##cputime=$(qstat -f -t "${SLURM_JOBID}" | grep 'resources_used.cput' | awk -F ' = ' '{print $2}')
##echo "#### cputime used by ${SLURM_JOBID} was: ${cputime:-null}" >> outputs/log.txt
####qstat -f -t ${SLURM_JOBID} >> outputs/log.txt


#!/bin/bash
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
###SBATCH -n 4 # number of cores
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G # memory pool for all cores
#SBATCH --time=24:00:00 # time (D-HH:MM)
#SBATCH --array=0-167

echo "In bash script, task = ${SLURM_ARRAY_TASK_ID}"

chr_list=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "X" "Y")

pseudo_ref_filepath="/projects/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8/"

founder_list=("A_J.39.fa" "C57BL_6J.39.fa" \
             "129S1_SvImJ.39.fa" "NOD_ShiLtJ.39.fa" \
             "NZO_HlLtJ.39.fa" "CAST_EiJ.39.fa" \
             "PWK_PhJ.39.fa" "WSB_EiJ.39.fa")

founder_index=$((${SLURM_ARRAY_TASK_ID}/21))
chr_index=$((${SLURM_ARRAY_TASK_ID}%21))
founder_ref=${founder_list[$founder_index]}
chr=${chr_list[$chr_index]}
echo $founder_index  $founder_ref  $chr_index  $chr

ref_fastscratch="/fastscratch/chenm/"$founder_ref

if [ ! -f $ref_fastscratch ]
then
    echo $ref_fastscratch" not found, copying $pseudo_ref_filepath to fastscratch"
    echo "cp $pesudo_ref_filepath$founder_ref $ref_fastscrach"
    cp $pesudo_ref_filepath$founder_ref $ref_fastscrach
fi

module load singularity
singularity exec  /projects/researchit/cervaf/containers/cnn_autoencoder_latest.sif python run_enformer.py $chr $ref_fastscratch

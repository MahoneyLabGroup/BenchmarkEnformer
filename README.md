# BenchmarkEnformer

run_enformer.py:   python script to generate Enformer predictions for genes on a given chromosome and a ref.fa file 

                   module load singularity
             
                   singularity exec  /projects/researchit/cervaf/containers/cnn_autoencoder_latest.sif python run_enformer.py $chr $ref_fa
                   
torun_enformer.sh: bash script to run run_enformer.py for 21 chromesomes (1-19,X,Y) of 8 founder pseudo references under:    
                    /projects/omics_share/mouse/GRCm39/genome/sequence/imputed/rel_2112_v8/                     
                    The transcription start site annotations are gtf files under:                    
                    /projects/compsci/omics_share/mouse/GRCm39/transcriptome/annotation/imputed/rel_2112_v8/                    
                    it will launch 168 (21 * 8) jobs 
                    
                    sbatch torun_enformer.sh

calc_qtl_coeff.Rmd: R script to calcualte 8 allel effects for each gene using the DO kidney dataset
                    output: data/coeff.csv

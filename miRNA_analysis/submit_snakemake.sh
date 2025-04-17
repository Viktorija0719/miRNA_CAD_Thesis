  #!/bin/bash
  #PBS -N snakemake_job
  #PBS -q batch
  #PBS -l walltime=48:00:00
  #PBS -l nodes=1:ppn=16,mem=160g
  #PBS -j oe
  #PBS -A rsu
  
  # Load the Conda module
  module load anaconda/conda-23.1.0
  
  # Activate the Conda environment where Snakemake is installed
  source activate snakemake_mirdeep2
  
  # Navigate to the directory containing the Snakemake workflow
  cd /home_beegfs/vikja01/miRNA_CAD_Thesis/miRNA_analysis
  
  # Execute Snakemake; adjust the number of cores and other parameters as needed
  
  snakemake -s Snakefile02 --cores 16 --use-conda
  
  # Deactivate the environment
  conda deactivate
  

# Sm_Mira_IvT
 For this project, our goal is to determine if miracidia from the liver and the intestines are transcriptomically different. If they are, we also want to identify any differentially expressed genes.

## Reduced Data vs Whole Data
 Since the data was very large, I did not have enough time to run it on the whole dataset. The reduced dataset has 10% of the data for each file. Anything with "_reduced" uses the reduced data.  
 The pipeline_star section is on the whole dataset, but only with the star alignment.

## Ways to Run the Analysis
 There are two ways of running this analysis: in a jupyter notebook or through snakemake.

### To Run Notebook:
 First create an environment with the env.yaml file. You will also need to run the line "BiocManager::install("GenomeInfoDbData", lib="/data/users/YOUR_USERNAME/.conda/envs/CONDA_ENVIRONMENT/lib/R/library/")" for the Rscript to work, changing the parts in all caps. Then run each part of the notebook.

### To Run Snakemake:
 First download snakemake. Then cd into pipeline or pipeline_reduced, depending on which you want to run. Finally, run snakemake with the --use-conda argument.
 
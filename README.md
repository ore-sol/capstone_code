## About the Pipeline
The purpose of this pipeline is to take paired end 16S data and generate genome scale metabolic models and associated matrices.  

## INPUT 
To start, a folder called "fastqfiles_1" is expected to be in the same directory as the 'pipeline_scripts' folder.  
Inside of the fastqfiles_1 folder are to be paired end 16S .fastq files. The file endings are expected to be '_R1_001' and 'R2_001'.  

## REQUIRED CONDA ENV
The script assumes a set of conda environments that contain the software called. These would be:  
1. DADA2  
```
conda create --name dada2-v1.40 bioconda::bioconductor-dada2
```  
2. NCBI  
```
conda create --name ncbi_download conda-forge::ncbi-datasets-cli
```  
3. Barrnap  
```
conda create --name barrnap bioconda::barrnap
```  
4. MASH  
```
conda create --name mash bioconda::mash
```  
5. Gapseq  
```
conda create --name gapseq bioconda::gapseq
```  
  
## After Creating Draft Models  
Given your samples, the script expects files in the format: SAMPLENAME_genera.csv that contains the genera identified in your sample (this can be done using phyloseq and then filtering)
        This will allow for separate creation of community interaction structure values
        Move these .csv files to the same folder that contains your model summary .csv files
Then:   
  
1. Run 'rename_genus_folders.sh' in same folder as the one that contains model summary .csv files (last thing generated in pipeline)
2. Run 'matrix_v6.R' in same folder as renamed models 
3. To calculate C, F and proxy-b values, run the following commands:   
    `source analyze_uptake_leakage_proxy_b_OS.R`  
    `run_batch()`  

## Additional Notes
This pipeline is run as sequential steps and users are invited to run either everything at once or individually to suit their project workflow. 



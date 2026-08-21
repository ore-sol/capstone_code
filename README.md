## About the Pipeline
The purpose of this pipeline is to take paired end 16S data and generate genome scale metabolic models and associated matrices.  

## INPUT 
To start, a folder called "fastqfiles_1" is expected to be in the same directory as the 'pipeline_scripts' folder.  
Inside of the fastqfiles_1 folder are to be paired end 16S .fastq files. The file endings are expected to be '_R1_001' and 'R2_001'.  

## REQUIRED CONDA ENV
The script assumes a set of conda environments that contain the software called. These would be:  
1. DADA2  
2. NCBI  
3. Barrnap  
4. MASH  
5. Gapseq  

## Additional Notes
This pipeline is run as sequential steps and users are invited to run either everything at once or individually to set their project workflow. 



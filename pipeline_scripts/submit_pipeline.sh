#!/bin/bash
# submit_pipeline.sh
JOB1=$(qsub step1_submit_dada2.pbs)
JOB2=$(qsub -W depend=afterok:$JOB1 submit_taxa.pbs)
JOB3=$(qsub -W depend=afterok:$JOB2 step2_accession_v5.pbs)
JOB4=$(qsub -W depend=afterok:$JOB3 step3_run_barrnap2.pbs)
JOB5=$(qsub -W depend=afterok:$JOB4 step4_extractandmash.pbs)
JOB6=$(qsub -W depend=afterok:$JOB5 step5_run_gapseq3.pbs)
JOB7=$(qsub -W depend=afterok:$JOB6 generate_summaries.pbs)
echo "Pipeline submitted. Final job: $JOB6"

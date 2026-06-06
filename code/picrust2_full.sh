
conda activate picrust2

# --- running the whole pipeline at once

#picrust2_pipeline.py -s ../Results/processed_sra_pipeline_v1/asv_sequences.fasta -i ../Results/processed_sra_pipeline_v1/asv_table.tsv -o picrust2_out_pipeline_split -p 1

place_seqs.py -s ../Results/processed_sra_pipeline_v1/asv_sequences.fasta -o bac_placed_seqs.tre -p 1 \
                    -r bacteria \
                    --intermediate placement_working_bac

hsp.py -i 16S -r bacteria -t bac_placed_seqs.tre -o bac_marker_nsti_predicted.tsv.gz -p 1 -n

hsp.py -i EC -r bacteria -t bac_reduced_placed_seqs.tre -o bac_EC_predicted.tsv.gz -p 1

hsp.py -i KO -r bacteria -t bac_reduced_placed_seqs.tre -o bac_KO_predicted.tsv.gz -p 1

combine_domains.py --table_dom1 bac_EC_predicted.tsv.gz --table_dom2 arc_EC_predicted.tsv.gz -o combined_EC_predicted.tsv.gz

combine_domains.py --table_dom1 bac_KO_predicted.tsv.gz --table_dom2 arc_KO_predicted.tsv.gz -o combined_KO_predicted.tsv.gz

metagenome_pipeline.py -i ../Results/processed_sra_pipeline_v1/asv_table.tsv \
                       -m combined_marker_nsti_predicted.tsv.gz \
                       -f combined_EC_predicted.tsv.gz \
                       -o EC_metagenome_out

metagenome_pipeline.py -i ../Results/processed_sra_pipeline_v1/asv_table.tsv \
                       -m combined_marker_nsti_predicted.tsv.gz \
                       -f combined_KO_predicted.tsv.gz \
                       -o KO_metagenome_out

pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o pathways_out \
                    --intermediate pathways_working \
                    -p 1

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
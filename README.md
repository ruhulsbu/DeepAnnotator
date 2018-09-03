# DeepAnnotator (Will be Updated on Sept 7, 2018 with Codebase and Dataset)

python general_lstm_gene_prediction.py start/stop codon

-->
time python general_lstm_start_classification.py start > train_start_codon.txt &

python evaluate_codon_classification_network.py model_name 0.50
python create_training_plot.py output_training_evaluation_testing_prediction.txt
python evaluate_codon_classification_network.py start general_lstm_model_start 0.50
python create_whole_genome_prediction.py gene_prediction_lstm_model 0.50 genome

-->

python create_whole_genome_prediction.py gene_prediction_lstm_model 0.50 genome
python create_whole_genome_prediction.py general_lstm_model_start 0.50 start
python create_whole_genome_prediction.py general_lstm_model_stop 0.50 stop

-->

python improved_whole_genome_prediction.py gzip_file_name
python intrag_coding_prediction_score.py gzip_file_name
python trying_nodp_solve_gene_boundary.py gzip_file_name

-->

python floatpt_whole_genome_prediction.py
python trying_dynp_solve_gene_boundary.py floatpt_file_name


python trying_dynp_solve_gene_boundary.py floatpt_whole_sequences_batch_4_3.gz > debug_dynp.txt




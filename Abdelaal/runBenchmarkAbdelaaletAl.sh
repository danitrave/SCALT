#! /bin/bash

#### Unzip files from original benchmark #### 
#unzip scRNAseq_Benchmark_datasets.zip

#### Group files in a directory ####
#mkdir ./originalFiles
#mv Inter-dataset/ Intra-dataset/ Rejection/ Scalability/ scRNAseq_Benchmark_datasets.zip originalFiles/

#### Transpose ####
#python3 transposition.py

#### Intra benchmark ####
#python3 intra_benchmark.py

#### Inter benchmark ####
python3 inter_benchmark.py

#### Classification with SCALT ####
#python3 classificationWithSCALT.py

#### Move all files in a directory ####
#mkdir ./AbdelaalEtAll_benchmark_02
#mv processedData/ originalFiles/ INTRA_benchmark/ INTER_benchmark/ intra_benchmark.py inter_benchmark.py medianF1_score.py cv_folds.R RdataParserCrossValidation.R RdataParser.R transposition.py classify_withSCALT/ AbdelaalEtAll_benchmark_02/


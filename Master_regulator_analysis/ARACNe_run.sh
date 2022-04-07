#Script - Running ARACNe using gene expression data
#Date of Creation - 01/24/2022

#Instructions
#Unzip the ARACNe.zip file and save a copy in the current directory
#Convert the gene expression data (raw counts) to tpm and save it as .rds object
#Run the following script
#Note: For samples with > 1000 cells, the script may require considerable amount of memory (10GB to 20GB) for the purpose of bootstrapping

#!/bin/bash

bash /Master_regulator_analysis/ARACNe/aracne-script.sh \
-n MBM05sn_ar \
-b /Master_regulator_analysis/MBM05_sn/ \
-i /Master_regulator_analysis/MBM05_sn.tpm.rds \
-r /Master_regulator_analysis/tf.cotf.sig.surf.SYMBOLS.txt

#After the script successfully runs, the following files would be created in the current directory:
# 1. MPM05sn_ar_ARACNe-table.tsv
# 2. MPM05sn_ar_finalNet-merged.tsv
# 3. MPM05sn_ar_unPruned.rds
# 4. MPM05sn_ar_pruned.rds
# 5. MPM05sn_ar_log-files
# 6. MPM05sn_ar_metaData.csv
# 7. MPM05sn_ar_net-tsvs

# For the purpose of running VIPER, one needs to select MPM05sn_ar_pruned.rds file (which is the pruned form of the interactome)

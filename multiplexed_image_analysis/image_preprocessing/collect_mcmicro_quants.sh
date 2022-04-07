## This script will copy all quantification files from their MCMICRO subfolders into one central folder for easy downstream analysis
project_dir="Melanoma_Brain_Metastasis"
main_mcmicro_folder="$project_dir/mcmicro"
quant_folder="$project_dir/mcmicro_quants"
for FILE in $main_mcmicro_folder/*
do
        if [ ! -f $quant_folder/$FILE ]
        then
                echo $FILE
                cp $FILE/quantification/* $quant_folder
        else
                echo $FILE" exists!"
        fi
done

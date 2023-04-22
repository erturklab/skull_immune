#! /bin/bash

#SBATCH -p cpu_p
#SBATCH -o /home/icb/%u/logs/cpdb_%j.out
#SBATCH -e /home/icb/%u/logs/cpdb_%j.err
#SBATCH -J a1_cpdb
#SBATCH -c 18
#SBATCH --mem=40G 
#SBATCH -t 16:00:00 
#SBATCH --nice=10000 

DATA_DIR='/home/louiskuemmerle/workspace/projects/A1/data2'
RESULTS_DIR='/home/louiskuemmerle/workspace/projects/A1/results'
DATASET='cellxgene_oct22_wSham_umaps'

echo 'start job'

source ~/miniconda3/bin/activate
source $HOME/.bashrc

echo 'activate environment a1_scanpy'
conda activate a1_scanpy
echo 'python script to save logits'
python ./script_cpdb_data.py $DATA_DIR $DATASET $1 $2 $3 $4 $5

echo 'deactivate environment a1_scanpy'
conda deactivate
echo 'activate environment a1_cpdb' # python 3.6 with cellphonedb 2.1.4
conda activate a1_cpdb


TEST=''
if [ $4 == 1 ]
then
    TEST='_test'
fi

echo $TEST

echo 'start cellphonedb analysis'

FOLDER_NAME=$1_$2_$3_seed$5$TEST

out_dir="${RESULTS_DIR}/cellphone_out_it1000/${FOLDER_NAME}/"
mkdir $out_dir
FILE1="${DATA_DIR}/cellphonedb/${FOLDER_NAME}/${DATASET}_non_log_meta.txt"
FILE2="${DATA_DIR}/cellphonedb/${FOLDER_NAME}/${DATASET}_non_log.txt"


echo "file1, file2, outdir:"
echo $FILE1
echo $FILE2
echo $out_dir

if [ -f "${out_dir}means.txt" ]; then
    echo "-- SKIP -- output in ${FOLDER_NAME} already exists."
else 
    cellphonedb method statistical_analysis $FILE1 $FILE2 --iterations=1000 --threads=16 --counts-data=gene_name --output-path=$out_dir
fi


#echo 'cellphone dotplot'
#
#cellphonedb plot dot_plot --means-path=$RESULTS_DIR/cellphone_out/$FOLDER_NAME/means.txt --pvalues-path=$RESULTS_DIR/cellphone_out/$FOLDER_NAME/pvalues.txt --output-path=$RESULTS_DIR/cellphone_out/$FOLDER_NAME/plots --output-name=dot_plot.pdf
#
#echo 'cellphone heatmap'
#mkdir $RESULTS_DIR/cellphone_out/$FOLDER_NAME/plots

#cellphonedb plot heatmap_plot $DATA_DIR/cellphonedb/$FOLDER_NAME/${DATASET}_non_log_meta.txt --pvalues-path=$RESULTS_DIR/cellphone_out/$FOLDER_NAME/pvalues.txt --output-path=$RESULTS_DIR/cellphone_out/$FOLDER_NAME/plots

echo 'finished'



CELLSNP_DIR="/storage/groups/ml01/workspace/louis.kuemmerle/projects/A1/results/donor_deconv/"
OUTPUT_DIR="/storage/groups/ml01/workspace/louis.kuemmerle/projects/A1/results/donor_deconv/vireo" 
JOBS_DIR="/home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/jobs_donor_deconv_vireo"

PARTITION="cpu_p"

mkdir $OUTPUT_DIR
mkdir $JOBS_DIR

SAMPLES="
    MUC12819
    MUC12820
    MUC12821
    MUC12822
    MUC12823
    MUC12824
    MUC12825
    MUC12826
    MUC12827
    MUC12828
    MUC12829
    MUC12830
    MUC12831
    MUC12832
    MUC12833
    MUC12834
    MUC12835
    MUC12836
    MUC12837
    MUC12838
    MUC12839
    MUC12840
    MUC12841
    MUC12842
    MUC12843
    MUC12844
    MUC12845
    MUC12846
    MUC12847
    MUC12848
    MUC12849
    MUC12850
"

for SAMPLE in $SAMPLES
do

sleep 0.1
job_file="${JOBS_DIR}/donor_deconv_${SAMPLE}.cmd"
cell_data="${CELLSNP_DIR}/${SAMPLE}"
out_dir="${OUTPUT_DIR}/${SAMPLE}"
    
echo "#!/bin/bash
#SBATCH -J donor_deconv_vireo_${SAMPLE}
#SBATCH -o ${JOBS_DIR}/donor_deconv_vireo_${SAMPLE}_%j.out
#SBATCH -e ${JOBS_DIR}/donor_deconv_vireo_${SAMPLE}_%j.out
#SBATCH -p ${PARTITION}
#SBATCH -t 3-00:00:00
#SBATCH -c 24
#SBATCH --mem=50G
#SBATCH --nice=10000

source ${HOME}/.bash_profile
conda activate a1_donor_deconv

mkdir $out_dir

vireo -c $cell_data -N 3 -o $out_dir

" > ${job_file}
sbatch $job_file

done

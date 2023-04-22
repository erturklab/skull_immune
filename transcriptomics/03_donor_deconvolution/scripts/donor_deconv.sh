

BAM_FILES_DIR="/storage/groups/bcf/projects/SingleCell/10x/projects/ertuerk/MCAO_vs_sham/rev5/cellranger"
OUTPUT_DIR="/storage/groups/ml01/workspace/louis.kuemmerle/projects/A1/results/donor_deconv" 
BARCODE_DIR="${OUTPUT_DIR}/barcodes"
JOBS_DIR="/home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/jobs_donor_deconv"

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

# Mouse has chromosomes 1-19 (+ X and Y, I think X,Y,MT are ignored) script fails if not defined 
# since it expects 22 chromosomes (human)
CHROMOSOMES="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y"   #"MT"  #"1,2"  #

for SAMPLE in $SAMPLES
do

sleep 0.1
job_file="${JOBS_DIR}/donor_deconv_${SAMPLE}.cmd"
out_dir="${OUTPUT_DIR}/${SAMPLE}"
bam="${BAM_FILES_DIR}/${SAMPLE}/bamfile/possorted_genome_bam.bam"
barcode="${BARCODE_DIR}/${SAMPLE}_no_QC.txt"
#barcode="${BARCODE_DIR}/${SAMPLE}.txt"
#barcode="${BAM_FILES_DIR}/${SAMPLE}/count_matrices/raw_feature_bc_matrix/barcodes.tsv"

    
echo "#!/bin/bash
#SBATCH -J donor_deconv_${SAMPLE}
#SBATCH -o ${JOBS_DIR}/donor_deconv_${SAMPLE}_%j.out
#SBATCH -e ${JOBS_DIR}/donor_deconv_${SAMPLE}_%j.out
#SBATCH -p ${PARTITION}
#SBATCH -t 3-00:00:00
#SBATCH -c 24
#SBATCH --mem=50G
#SBATCH --nice=10000

source ${HOME}/.bash_profile
conda activate a1_donor_deconv

mkdir $out_dir

cellsnp-lite -s $bam -b $barcode -O $out_dir -p 22 --chrom $CHROMOSOMES --minMAF 0.1 --minCOUNT 50 --gzip
" > ${job_file}
sbatch $job_file

done

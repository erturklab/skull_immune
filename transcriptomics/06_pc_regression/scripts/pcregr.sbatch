#!/bin/bash

#SBATCH -o /home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/pcregression/logs/out_%j.txt
#SBATCH -e /home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/pcregression/logs/out_%j.txt
#SBATCH -J PCR
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 0-02:00:00
#SBATCH --nice=10000
#SBATCH -p cpu_p

source $HOME/.bashrc

conda activate a1_img

WD="/home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/pcregression/"

python "${WD}pcregr.py" $1 $2


#!/bin/bash

WD="/home/icb/louis.kuemmerle/projects/a1/code/skull_immune/scripts/pcregression/"

declare -a CTS=(
    "progenitors"
    "neutrophil"
    "monocyte"
    "B cell"
    "T cell"
    "NK cell"
    "NK-T cell"
    "dendritic cell"
    "macrophage"
    "microglia"
    "erythroid precursor"
    "erythroid cell"
    "basophil"
    "structural cell"
    "brain cell"
    "megakaryocyte"
    "innate lymphoid cell"
)


for CT in "${CTS[@]}"
do
    CT=${CT// /_} # replace all spaces with underscores (in pcregr.py they are converted back)
    sbatch "${WD}pcregr.sbatch" level1 $CT
    sleep 1

done



#declare -a CTS=(
#    "hematopoietic stem cell"
#    "common myeloid progenitor"
#    "granulocyte-monocyte progenitor"
#    "neutrophil-primed GMP"
#    "monocyte-primed GMP"
#    "erythroid progenitor"
#    "common DC progenitor (CDP)"
#    "pro neutrophil"
#    "pre neutrophil"
#    "immature neutrophil"
#    "mature neutrophil"
#    "monocyte progenitor"
#    "classical monocyte"
#    "non-classical monocyte"
#    "pro B cell"
#    "pre B cell"
#    "immature B cell"
#    "mature B cell"
#    "plasma cell"
#    "Cd8 T cell"
#    "Cd4 T cell"
#    "gdT cell"
#    "NK-T cell"
#    "NK cell"
#    "plasmacytoid DC"
#    "conventional DC1"
#    "conventional DC2"
#    "perivascular macrophage"
#    "monocyte-derived macrophage"
#    "microglia"
#    "macrophage"
#    "antigen-presenting macrophage"
#    "erythroid cell"
#    "erythroblast"
#    "erythrocyte"
#    "basophil progenitor"
#    "basophil"
#    "fibroblast"
#    "dural fibroblast"
#    "endothelial cell"
#    "adipose-derived stromal cell"
#    "brain-Chroid Plexus endothelial cell"
#    "meningeal-Choroid Plexus cell"
#    "Omp+ cell"
#    "Gnb3+ cell"
#    "neuron"
#    "astrocyte"
#    "oligodendrocyte"
#    "megakaryocyte"
#    "innate lymphoid cell"
#)
#
#
#for CT in "${CTS[@]}"
#do
#    CT=${CT// /_} # replace all spaces with underscores (in pcregr.py they are converted back)
#    sbatch "${WD}pcregr.sbatch" level2 $CT
#    sleep 1
#
#done
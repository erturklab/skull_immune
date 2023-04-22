#! /bin/bash

TEST=0
SEEDS=(0 1 2 3 4)
# seed 0 is with gemoetric sketching, any other seed is with random sampling --> was done to check for robustness of results

for CT_SELECTION_KEY in level1 level2_min10
do
    for SEED in ${SEEDS[@]}
    do
    	for CONDITION in Naive Sham MCAO
    	do
       	    for BONE in Skull Pelvis Femur Vertebra Humerus Scapula
        	do
                #sbatch script_cpdb.sh $CT_SELECTION_KEY $CONDITION $BONE $TEST $SEED
            	bash script_cpdb.sh $CT_SELECTION_KEY $CONDITION $BONE $TEST $SEED
            	sleep 2
            done
        done
    done
done

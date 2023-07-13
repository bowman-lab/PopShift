#!/bin/bash

#calx=~/bowmore/software/MSM-docking/msm_binder.py
calx=~/python_modules/MSM-docking/msm_binder.py
#extract=~/bowmore/software/MSM-docking/extract-scores.py
extract=~/python_modules/MSM-docking/extract-scores.py
ncores=12
for runname in {15..25..5}xsmina; do
    for fpd in {dt-100-100,dt-400-50,rmsd-50,rmsd-200}-l-10-p-3; do
            rundir=$fpd
            dock_run=$rundir/$runname
            extracted_scores=$rundir/extracted_scores/$runname
            eq_probs=$fpd-eq_probs.npy
            # python extract-eq-emma.py $model.pickle $model-eq-probs.npy
            python $extract -n $ncores -t smina $dock_run
            ls $extracted_scores/*.h5
            python $calx\
             -n $ncores\
             --reweighted-eq\
             --out $extracted_scores/binding-calx\
             bin-samples\
             $eq_probs\
             $extracted_scores/*.h5
    done
done

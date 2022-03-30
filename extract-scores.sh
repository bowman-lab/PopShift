#!/bin/bash

extract_scores() {
    reference_path=$1
    ref_mol_name=$2
    suffix=$4
    # clear out and write outfiles
    for trj in $reference_path/*; do
        rep_counter=0
        for rep in $trj/$ref_mol_name/*; do
            outfn=$3-$rep_counter.dat
            ordered_path_list=$3-$rep_counter.out
            if [[ -f $ordered_path_list ]]; then
                rm $ordered_path_list
            fi
            echo '# trj frame score' > $outfn
            rep_counter=$((rep_counter + 1))
        done
    done

    trj_counter=0
    for trj in $reference_path/*; do
        rep_counter=0
        outfn=$3-$rep_counter.dat
        ordered_path_list=$3-$rep_counter.out
        echo $trj
        echo $trj/$ref_mol_name
        for rep in $trj/$ref_mol_name/*; do
            frame_counter=0
            echo $rep
            for frame in $rep/frame*$suffix; do
                score=$(awk '/REMARK VINA RESULT/ {print $4; exit}' $frame)
                echo $trj_counter $frame_counter $score >> $outfn
                echo $frame >> $ordered_path_list
                frame_counter=$((frame_counter + 1))
            done
            rep_counter=$((rep_counter + 1))
        done
        trj_counter=$((trj_counter + 1))
    done
}
export -f extract_scores
refpath=myh2
suffix=44.65_91.42_53.83_14.0_12.0_18.0
scoresuff=$suffix
framesuff=$suffix.pdbqt
#for mol in blebbistatin mt-{100..105} mt-{111..115} mt-{131..135}; do
##for mol in blebbistatin; do
#    extract_scores $refpath $mol $mol-scores
#done
for suffix in 44.65_91.42_53.83_14.0_12.0_18.0 38.65_93.42_53.83_24.0_18.0_22.0; do
    scoresuff=$suffix
    framesuff=$suffix.pdbqt
    parallel --jobs 75% extract_scores $refpath {} {}-scores-$scoresuff $framesuff ::: \
        blebbistatin mt-{100..105} mt-{111..115} mt-{131..135}
done

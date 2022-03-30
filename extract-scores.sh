#!/bin/bash


# This function is hard-coded to interact with the directory structure created by the docking parallel tools.
# While its structure might be useful inspiration if you reorganized your data directories, it is unlikely to work
# in the context of any real changes.
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
            # create the outfile by (over)writing a header line
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
# first argument to this script is the reference msm path
refpath=$1
# shift indexing of arguments by one, to reflect that we've dealt with arg 1 already
shift 1
suffix=$1
shift 1
scoresuff=$suffix
framesuff=$suffix.pdbqt
# use remaining args as compound names, run parallel extraction jobs for each
parallel --jobs 75% extract_scores $refpath {} {}-scores-$scoresuff $framesuff ::: $@

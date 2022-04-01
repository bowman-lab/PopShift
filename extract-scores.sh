#!/bin/bash


# This function is hard-coded to interact with the directory structure created by the docking parallel tools.
# While its structure might be useful inspiration if you reorganized your data directories, it is unlikely to work
# in the context of any real changes.
extract_scores() {
    reference_path=$1
    ref_mol_name=$2
    outfile_prefix=$3
    suffix=$4
    awk_script=$5

    # clear out and write outfiles
    for trj in $reference_path/*; do
        rep_counter=0
        for rep in $trj/$ref_mol_name/*; do
            outfn=$outfile_prefix-$rep_counter.dat
            ordered_path_list=$outfile_prefix-$rep_counter.out
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
        outfn=$outfile_prefix-$rep_counter.dat
        ordered_path_list=$outfile_prefix-$rep_counter.out
        echo $trj
        echo $trj/$ref_mol_name
        for rep in $trj/$ref_mol_name/*; do
            frame_counter=0
            echo $rep
            for frame in $rep/*$suffix; do
                # use memory-mapped awk for better performance
                score=$(mawk -f $awk_script $frame)
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
# if this string is provided to the program as several suffixes
# within quotes and separated by whitespace,
# the script will iterate over each
suffixes=$1
shift 1
awk_script=$1
# use heredocs to write awkscripts, so that there's no trouble with bash preemptively expanding awk variables.
if [[ $awk_script == 'vina' ]]; then
cat << 'EOF' > $awk_script.awk
/REMARK VINA RESULT/ {print $4; exit}
EOF
awk_script=$awk_script.awk
elif [[ $awk_script == 'smina' ]]; then
cat << 'EOF' > $awk_script.awk
/Result/ {print $3; exit}
EOF
awk_script=$awk_script.awk
else
    echo "Third argument (echoed below) is not one of 'smina' or 'vina'; interpreting as an awk command file name."
    echo "$awk_script"
    awk_script="$awk_script"
fi
shift 1
for suffix in $suffixes; do
    scoresuff=$suffix
    framesuff=$suffix.pdbqt
    # use remaining args as compound names, run parallel extraction jobs for each
    parallel --jobs 100% extract_scores $refpath {} {}-scores-$scoresuff $framesuff "$awk_script" ::: $@
done

# introduction

this is a repository containing some command-line tools for implementing the popshift framework for docking. the purpose behind these tools are to allow users both to make alterations to these workflows for their own use, and also to apply them in similar fashion to the use they were put toward in the original [popshift manuscript](https://www.biorxiv.org/content/10.1101/2023.07.14.549110v2). they're provided with the gplv2 license, in the hope that others find them useful in their own projects, but also with the explicit intention of receiving improvements as others make them in their own contexts. below you'll find a discussion of the general workflow from that paper, as well as some documentation of what each script does. note that these command line tools all also have `argparse` i/o tools, so if you're wondering about how to use them in more detail calling them with no arguments, or with the `-h` or `--help` flags will print some help text that may be more rapidly useful than leafing through this document.

# Installation

these python scripts are designed to be used as command-line tools; some can also be imported into python interpreter sessions or notebooks, if one is interested in calling the functions from them independently. thus, getting python environments set up with the correct modules present is the main installation task. getting the parameterization program associated with the docking code you're interested in installed is the other major task.

### Minimal install using `smina` for docking and Openbabel for system prep.
```
# create the analysis postprocessing and framepicking environment
mamba create -n popshift -c conda-forge obabel loos cython mdtraj 
mamba activate popshift
git clone https://github.com/bowman-lab/enspara
cd enspara
python setup.py install
# create the docking environment
mamba create -n smina -c conda-forge smina obabel
```

Because of artifical conflicts in `boost` version `loos` and `smina` are currently mutually exclusive in the same env. If either were built as compiled software this issue would not exist, as both are in fact compatible with a wide range of `boost` versions.

## System preparation programs

the two preparation programs we've used are for autodock family programs, and generate pdbqt files. they are:

### Openbabel

get it from conda-forge:
`mamba install -c conda-forge openbabel`

get it from apt (for debian-like systems):
`sudo apt install openbabel`

### `prepare_receptor` and `prepare_ligand`:

because these programs _require_ python 2, and therefore don't play nicely with any other python tools written anytime recently, they are best installed from the tarball available from the [adfr website](https://ccsb.scripps.edu/adfr/downloads/). you should follow the directions there to get them added to your path.

note that you can view the [licenses](https://ccsb.scripps.edu/adfr/license/) these programs are being provided under on that same site. because we only need `prepare_receptor` and `prepare_ligand`, these are available under lgplv2 at time of writing.

## Docking programs

right now our scripts are set up to either use [vina](https://vina.scripps.edu/) or [smina](https://github.com/mwojcikowski/smina). in both cases the easiest way to get either set up is to install into a conda env for docking:
```
mamba create -n docking -c conda-forge smina vina
```

we have run into an issue with using vina from `conda-forge` source; in particular, we found [there was a bug](https://github.com/ccsb-scripps/autodock-vina/issues/90) that produced pretty significant variability in docking scores given the same pose, which is pretty bad for this application where the scores themselves are being used. we also noticed another bug whereby sometimes, seemingly at random across thousands of docking runs, a handful of poses will be given an extra 10-15 kcal/mol of favorability. there are several other energy function bugs on their issue tracker that this could be. until the version on conda forge is greater than `1.2.4` we'd recommend compiling the `develop` branch of [vina hosted on github](https://github.com/ccsb-scripps/autodock-vina), following their install instructions for [building from source](https://autodock-vina.readthedocs.io/en/latest/installation.html#building-from-source), focusing on the python bindings since that's the only part of the install used by `docking_parallel.py`.

note if you are trying to get the right compilers available you can install those from conda forge via the `compilers` repository. here's an example install line:
```
mamba create -n built-vina -c conda-forge compilers boost-cpp swig numpy
```
then you'd set up the python bindings:
```
conda activate built-vina
cd path/to/cloned/autodock-vina/build/python
rm -rf build dist *.egg-info
python setup.py build install
```

## Trajectory manipulation and postprocessing

to run the python tools associated with trajectory manipulation and post processing you'll need [`loos`](https://github.com/grossfieldlab/loos), a trajectory handling library, and [`enspara`](https://github.com/bowman-lab/enspara), an msm/utilities library. you'll also very likely need either [`pyemma`](http://www.emma-project.org/latest/), [`deeptime`](http://www.emma-project.org/latest/), or both.

a tidy way to install all of these into the same environment is to grab the libraries that are available on `conda-forge` and also the depeendencies that are needed for enspara, then clone enspara and build it using its `setup.py` script.

```
mamba create -n popshift -c conda-forge loos deeptime pyemma mdtraj cython mpi4py
conda activate popshift
git clone https://github.com/bowman-lab/enspara
cd enspara
python setup.py install
```

## Dependencies summary
1. ADFR suite (for `prep_parallel.py` and `prepare_ligand.py`)
    - > https://ccsb.scripps.edu/adfr/downloads/  
    cp prepare_ligand and prepare_receptor into $PATH
2. Autodock Vina python bindings (for Vina docking)
    - > git clone https://github.com/ccsb-scripps/AutoDock-Vina  
    cd AutoDock-Vina/build/python  
    conda install -c conda-forge numpy boost-cpp swig  
    rm -rf build dist *.egg-info (to clean previous installation)  
    python setup.py build install  
3. Smina
    - > mamba install -c conda-forge smina
4. enspara (for `popshift.py`, `pick_align_frames.py` `extract_scores.py`)
   - Follow [enspara](https://enspara.readthedocs.io/en/latest/installation.html) install instructions with `conda`.
   - If using PyEMMA discretized trajectories, then it must be installed:
     - > mamba install -c conda-forge pyemma
5. LOOS (for `pick_align_frames.py`, `draw_box.py`, `add_bonds_posthoc.py`, and `rmsd_receptor_ligand.py`)
    - > mamba install -c conda-forge loos
    - Can also be compiled following the instructions on the [loos github](https://github.com/GrossfieldLab/loos)

# Workflows

the scripts provided here are supposed to be modular enough that you can adjust them to your needs, or swap one part out in place of another part--we are certainly going to do this as our research on this framework evolves, so why shouldn't you? however, there are a few pathways that probably make general sense to understand, and the easiest way to begin to see what order things need to happen in is probably with some prescribed workflows. find some below:

## Docking to an msm

in the original popshift ms, we considered how to apply the popshift framework to the ensemble docking problem. the general workflow to actuate that study, with the tool(s) needed, follows:

1. obtain a satisfactory msm representing ligand-free simulations of your receptor.
2. determine a selection that picks out your pocket of interest.
3. pick and align frames from the msm using `pick_align_frames.py` to ensure that all the residues that should be in your box are nicely co-aligned.
   - check that things landed where you expected them to land using a visualizer like `pymol` and `draw_box.py`.
   - noting that the alignment will recenter the coordinate system of all the models at the centroid of the atoms selected as the 'pocket' atoms, so boxes centered at or close to `0,0,0` are probably good choices in general.
4. prepare receptors and ligands (using either openbabel or autodock tools `prepare_receptor` and `prepare_ligand`), here you have some alternatives:
   - use `prepare_ligand` and `prepare_receptor` on the command line directly with your favorite shell concurrency tool such as `xargs` or gnu `parallel`.
   - alternatively, use `obabel` with the `-o pdbqt` flag.
5. dock to each sample, saving the docking score in the file containing the ligand pose, using `docking_parallel.py`.
6. extract scores from these output files and collate them into arrays using `extract_scores.py`
7. use the score arrays to do popshift calculations, with `popshift.py`
8. plot stuff, compare to experiments, and generally *go nuts*, by reading the json that `popshift.py` outputs into your downstream scripts or notebooks.

note that the docking here is done with either vina or smina, so steps 4, 5, and 6 are somewhat particular to how those codes do that. you can of course color outside of the lines here, but you'll probably be adding some functionality to the scripts we have in order to do so.

### frame picking using `pick_align_frames.py`

for the docking below to work as completely independent jobs, we extract frames representing each state, then save them into subdirectories where each subdirectory name corresponds to the msm state index. each file name consists of a trajectory index, then a frame index, separated by a dash. here is an example of frames picked from a 25 state msm built from five longish trajectories (obtained by calling `tree -p '*.pdb' -v receptor`) output abbreviated for clarity:

```
receptor
├── 0
│   ├── 0-55257.pdb
│   ├── 0-56428.pdb
│   └── 1-43334.pdb
├── 1
│   ├── 2-105110.pdb
│   ├── 2-175513.pdb
│   └── 3-48068.pdb
├── 2
│   ├── 0-157652.pdb
│   ├── 0-158908.pdb
│   └── 0-191040.pdb
...
└── 24
    ├── 4-24060.pdb
    ├── 4-186351.pdb
    └── 4-189200.pdb

25 directories, 75 files
```
the commandline i used to create this was:

```
frame_picker=path/to/popshift/pick_align_frames.py

python $frame_picker \
  --find-bonds 1.95,1.35\
  --make-receptor-sel-chain-a\
  --assignments dtrajs.h5\
  --align-resid-list pocket-resids.txt\
  --number-frames 3\
  --mapping mapping.json\
  $system_name\
  prot_masses.pdb\
  eq-probs.npy\
  random\
  '!hydrogen'\
  traj_list.txt
```

here are the options and arguments, in order:
- `find-bonds`: find bonds using two different cutoffs, the first for heavy atoms and the second for hydrogen. not necessary if bonds are present in the file used for 'model'. better to have bonds written into pdbs so that the prep program doesn't need to build them.
- `make-receptor-sel-chain-a`: make the written pdbs have chain a identification. some system prep tools require this.
- `assignments`: the assignments or discretized trajectory file, saved as a 'ragged array. alternatively, an array of center indices saved as a `.npy`. to subsample randomly the full assignments file is obligatory, and **the order of the trajectories in the assignments file and the coordinate trajectory paths must match**.[^1]
- `align-resid-list`: this expects a text-file that `numpy.genfromtxt` can turn into a 1-d array of `int`s. these correspond to the `resid`s from the model that should be used in the iterative alignment. 
- `number-frames`: select this many frames from each msm bin.
- `mapping`: this corresponds to an enspara-style 'to' trim mapping; a dictionary containing the indexes of the original assignments trajectories as keys and the indexes they should be mapped to after ergodic trimming as values.
- `system_name`: the name you'd like to give to this particular extraction. will become the directory name within which the `receptor` dir containing the tree of selected frames is written.
- `prot_masses.pdb`: this is the model file. for this particular command i was using a `prot_masses.pdb` generated by gromacs, which means it doesn't have connectivity which is why the example also contains `find-bonds`.
- `eq-probs.npy`: a numpy array containing equilibrium probabilities. could also be a pyemma msm saved as a pickle file, in which case `mapping` would not have been necessary. used to communicate to the script how many states are active.
- `random`: the sampling mode requested here randomly picks frames from each msm bin without replacement. there are other options; consult the help statement. `centers` will switch the tool to operating on the provided coordinate trajectory as a centers trajectory.
- `!hydrogen`: the [loos selection string](https://grossfieldlab.github.io/loos/selections.html) to apply to the residues selected as the pocket for the purposes of alignment. this string means 'not hydrogen', or in other words 'all heavy atoms'.
- `traj_list.txt`: a text file containing paths to the coordinate trajectoies frames should be extracted from. a sequence of at least one paths may also be provided here, but as previously mentioned **the order of the trajectories in the assignments file and the coordinate trajectory paths must match**. because of this i much prefer to create a list of trajectories and save it as a text file when i cluster, then pass that list around to subsequent tools. if operating in centers mode then this argument should be a path to the trajectory of centers.

calling `pick_align_frames.py --help` will provide more information about alternative strategies for picking frames (for example, using just cluster centers, or picking particular numbers of frames from each bin). it's also worth noting that pick_align_frames can be imported into other python sessions, so if you want to rip random frames into sub-trajectories for other reasons this is possible. finally, note that subsample trajectories that are not pdbs can also be written to disk if desired. this can be nice if trying to use normal trajectory analysis tools on the selected frames, because reading in piles of pdbs instead of one dcd per msm bin can be annoying and slow.

[^1]: the order of `*.xtc`, python's `glob.glob('*.xtc')` python's `pathlib.path.glob('*.xtc')` and `find . -name '*.xtc'` will not necessarily be consistent across runs. the first of these results in ascii-betical ordering (meaning that `traj-11.xtc` will come _before_ `traj-2.xtc` for bash and derivative shells). the other three, if unsorted, will come in the file-system's order. importantly, **the filesystem order can change** if the contents of the filesystem change, which is why many references point out that the returns of such commands--when unsorted explicitly--are arbitrary. this is an easy way for one to shoot oneself in the foot.
# Docking scripts

## Preparing receptors and ligands for docking

we prepare each of the receptor conformations as a completely independent system, so each is represented by its own pdbqt file and can thus be run as an isolated job on whatever hardware is available. this leaves some things to be desired, but as a first pass it's worked for us.

in general these scripts are set up to play nicely with one another by assuming that the directory structure of the intermediate files will be maintained--this isn't a hard requirement, but it does keep things convenient. because we are using other command-line tools to prepare ligands and receptors for docking, we are really just expecting that the prepared files will be in the same place as the files they were prepped from. there are a lot of ways to call a command-line tool on a file such that it writes a new file with the same path but a different file extension. here are several i have used in the past:

using gnu parallel:

```
# use find to get the paths rel. to cwd;
# write to file to prevent arg-list too long
find frame_picking_system_dir/receptor/ -name '*.pdb' > receptor_paths.txt
# don't use too many jobs; i/o intensive task will be bound by disk perf.
parallel -j 8 obabel -opdbqt {} -o {}qt :::: receptor_paths.txt
```

or with `prepare_receptor`:

```
parallel -j 8 prepare_receptor -r {} -o {}qt :::: receptor_paths.txt
```

what's going on here? scrutinizing the docs for [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) suggests that parallel will be running the commandline on each of the elements in the file when in `::::` input mode. it substitutes each line into the 'command line' you provide as an argument wherever the curly braces are, then runs that as its own subprocess.

using [xargs](https://manpages.org/xargs) as an alternative to parallel:

```
find frame_picking_system_dir/receptor/ -name '*.pdb' | xargs -i% -n 1 -p 8 obabel -opdbqt % -o %qt
```

if only dealing with a few frame-picking runs (or maybe just one) for a relatively parsimonious msm (one with few states), we can afford to do this sequentially (nicer if you're getting errors, as well). in this case, by far the easiest thing to do is write a simple loop in the shell of your choosing. here's a while loop that takes advantage of the redirect-to-file from `find` given in the `parallel` example:

```
while read receptor_path; do
    echo $receptor_path
    obabel -opdbqt $receptor_path -o ${receptor_path}qt
done < receptor_paths.txt
```

### preparing ligands

ligands can also be prepped by command-line tools; the docking script `docking_parallel.py` expects you to provide the ligands as either a text-file of ligand paths or as a series of command-line arguments. as such, you can organize these files how you see fit. i usually put them in their own directory, then prep 'everything' with a certain file extension within that directory. so for example:

```
# assume there's a directory called ligands with .sdf or .mol2 files in it
parallel -j 8 obabel {} -h -isdf -opdbqt -o{.}.pdbqt ::: ligands/*.sdf
```
note: the `-h` flag adds hydrogens. don't throw it if you don't want that.

alternatively, if you have the adfr suite installed and would like to use `prepare_ligand` instead of openbabel (note you need `.mol2` or `.pdb` files for `prepare_ligand`, as of this writing):

```
parallel -j 8 prepare_ligand -l {} -ahydrogens -o{.}.pdbqt ::: ligands/*.mol2
```

## docking using jug

again here there are many ways to organize docking so long as the docking code is available on the command-line. we like to use jug, because it integrates well with our workload manager (slurm). we have also used it with lsf. following is an example sbatch script, with comments about the various options for docking paprallel.
```
#!/bin/bash
#sbatch --array=0-99
#sbatch -n 1
#sbatch -p all
#sbatch -o dockouts/%a-%a.out
#sbatch -e dockouts/%a-%a.out

receptor_dir=system_name_from_framepicking
jug execute -- path/to/popshift/docking_parallel.py\
 -d smina\
 -e 32\
 $receptor_dir/receptor\
 $receptor_dir/'12xsmina'\
 '0,0,0'\
 '12,12,12'\
 all-ligands/*.pdbqt

```
the sbatch pragma lines are just asking for a 100 element array job where each job is a one-core job with default everything on the 'all' partition. obviously you can adjust that if you'd like, but the tool expects to just run single-threaded docking jobs, so there's no point to increasing the per-job core count with the script as written. the command lines are as follows:
- `jug execute --`: saying this invokes jug instead of the python interpreter, issuing the execute command. the `--` tells the jug program that the part of the command line it should read is concluded; if this is omitted and the command line contains long-option flags that should be going to `docking_parallel.py` jug will interpret them as being for it instead.
- `d` is the 'docking-algorithm' flag, and picks from the available docking options--at present either `smina` or `vina`.
- `e` is the 'exhaustiveness' flag; for vina, the exhaustiveness flag is the number of mc runs to perform.
- `$receptor_dir/receptor`: the directory containing the receptor conformations to dock to.
- `$receptor_dir/12xsmina`: the nickname to give the docking run. it's preferable to keep all the docking runs organized in the same directory to their `receptor` dir for ease of later scripting. 
- `0,0,0`: the coordinates of the center, separated by commas. note that if you used `pick_align_frames.py` to get receptor structures to dock to, 0,0,0 is a relatively good choice for a center because the coordinates get zeroed by the frame 
- `12,12,12`: the edge lengths of the x, y, and z edges of the docking box.
- `all-ligands/*.pdbqt`: the last n arguments should be prepped ligand pdbqts. i bunged the prepped pdbqts of the ligands into a directory called `all_ligands`, then provided them to the command-line as a glob. note that which order the ligands get docked in is not in general important.

### output

the output from docking will be written into a directory structure that is within a directory given to the nickname for the docking run, with sub-directories using the stem of the file-names of the ligand pdbqts, and the directory structure beneath that containing numbered directories for each msm bin, and each docked pose containing _only_ the ligand coordinates will have a file-name matching the receptor structure to which it was docked. so for the above frame-picking example, supposing the ligand files docked in were `benzene.pdbqt`, `toluene.pdbqt`, and `p-xylene`:

```
tree -d 12xsmina
├── benzene
├── p-xylene
└── toluene

3 directories
tree 12xsmina/benzene
benzene
├── 0
│   ├── 0-55257.pdbqt
│   ├── 0-56428.pdbqt
│   └── 1-43334.pdbqt
├── 1
│   ├── 2-105110.pdbqt
│   ├── 2-175513.pdbqt
│   └── 3-48068.pdbqt
├── 2
│   ├── 0-157652.pdbqt
│   ├── 0-158908.pdbqt
│   └── 0-191040.pdbqt
...
└── 24
    ├── 4-24060.pdbqt
    ├── 4-186351.pdbqt
    └── 4-189200.pdbqt

25 directories, 75 files
```

 ### issuing jug commands
there are other [jug commands](https://jug.readthedocs.io/en/latest/subcommands.html). in order to issue them for a jug file that has command line arguments, one must make sure that the same command-line arguments are provided. one must also be using the same versions of python and jug, so its important to have the environment the dock job was launched from active. thus, often the easiest thing to do is make the sbatch into a script that takes an argument, and inserts that argument in place of `execute`:

```
#!/bin/bash

receptor_dir=system_name_from_framepicking
jug $1 -- path/to/popshift/docking_parallel.py\
 -d smina\
 -e 32\
 $receptor_dir/receptor\
 $receptor_dir/'12xsmina'\
 '0,0,0'\
 '12,12,12'\
 all-ligands/*.pdbqt
```
call this script `issue-jug-cmds-dp.sh` or similar. i could then use this script to obtain jug status reports, or to invalidate the results of a calculation so that i could rerun it:
```
issue-jug-cmds-dp.sh status
issue-jug-cmds-dp.sh 'invalidate dock_smina'
issue-jug-cmds-dp.sh cleanup
```

## analyze results with popshift framework

now we are ready to do something with the per-state affinity and pose estimates we've made. we'll extract the scores, then use them to make affinity estimates.

### extract scores

extract scores with `extract_scores.py`. this will pull the values of docking scores out of the various pose files and save them to an [`enspara`](https://enspara.readthedocs.io/en/latest/installation.html) ragged array in h5 format, since it is not in general true that the number of frames extracted per bin is always the same. this is a pretty simple tool; here is an example command line:

```
python path/to/popshift/extract_scores.py -n 8 -t smina 12xsmina
```
the commandline works as follows:
- `n`: the number of cores to use; this job is pretty i/o bound, but i have found that using a handful of cores can speed it up. 
- `t`: the type of output file the script can read. will match the `d` option from `docking_parallel.py`.
- `12xsmina`: the docking run nickname we would like to run the extraction out of.
  
by default, the resulting score extracts will be written to `extracted_scores`, in the same directory as the docking run is in. it takes the path to that directory, then goes up one level. within `extracted_scores` there will be a dir for each docking run you've extracted, and it will be named to match the directory name for the docking run. there will be an .h5 file within this directory for each ligand name. the elements will correspond to the scores extracted from the pose files.

### calculate macroscopic binding constants and reweighted state probabilities

to obtain binding free energy estimates and reweighted state populations, use `popshift.py`. this program was designed to either be a commandline tool that will perform these analyses on files arranged in the expected way, or as a small module with functions to perform the desired operations, if working within a notebook is desired. here's an example of command-line usage, but as with all the other tools this one has broader summary documentation provided by the output of the `--help` flag.

```
extracts=extracted_scores/12xsmina
python path/to/popshift/popshift.py \
  -n 8\
  --out $extracts/binding-calx\
  --reweighted-eq 1\
  bin-samples\
  path/to/eq-probs.npy\
  $extracts/*.h5
```

the commandline works as follows:
- `n`: number of threads--this code has `python mutltiprocessing` style threads, but unless one has huge numbers of ligands tasks are so fast using this threading is imperceptably different.
- `out`: the output path to write the affinity calculation json and any reweighted free energy (ragged) arrays to.
- `reweighted-eq`: the concentration of ligand you'd like to reweight your equilibrium probabilities by; here i've chosen 1m to simulate 'saturating' ligand conditions. note that reweighted probabilities will be per frame sampled, not per bin--they'll match the shape of the extracted scores arrays.
- `bin-samples`: the mode frame picking was performed in. should match the argument to `pick_align_frames.py`.
- `path/to/eq-probs.npy`: the equilibrium probabilities for the apo simulations used throughout.
- `$extracts/*.h5`: the extracted scores, in no particular order. assumes the stem of the path to each score `.h5` is the name of the ligand, for keying and naming its output.

Because the output is matched to the extracted scores, I like to write to a sub-directory of the extracted scores directory. If reweighted free energies were not requested, then only the scores post processing will be written. Three types of files may be created inside the directory provided as an arg to `--out`. One `calx.json` with the results of the affinity analysis inside it; one `{ligand_name}-eq_probs.h5` file with reweighted equilibrium probabilities corresponding to each _conformation_ analyzed. This array will have the same shape as the extracted scores arrays. There will also be a difference in population array, saved as a change in free energy per state between the frame weights and the reweights. In other words a $\Delta G$ per state from the initial population to the reweighted population. These values will be saved in a file called `{ligand_name}-dg.h5`. 

The output of the script will be a json with several different averages reported for each ligand, as well as a `"log"`. it's been indented to make it more human readable, but the standard python json module can turn it into a dictionary of dictionaries and lists with no special commands, so you can also inspect or plot its contents in a python interpreter or notebook as you see fit.  

Within the `"results"` object, there will be ligand names associated to sub-objects, each of which will be a dictionary keyed by ligand names with values that are themselves dictionaries with each of the type of calculations pop-shift can use to aggregate a binding site into a macroscopic dissocation constant. These are the `"popshift dG"` free energy with the provided kT, a `"popshift K_D"` which is the same popshift estimate of binding affinity as a dissociation constant in units of micromolar, the `"best score"`, the `"simple average"`, and the `"weighted average"`, which corresponds to the expectation over apo equilibrium probabilities we previously referred to as [boltzmann docking](https://www.nature.com/articles/ncomms12965).



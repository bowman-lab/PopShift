# Introduction

This is a repository containing some command-line tools for implementing the PopShift framework for docking. The purpose behind these tools are to allow users both to make alterations to these workflows for their own use, and also to apply them in similar fashion to the use they were put toward in the original [PopShift manuscript](https://www.biorxiv.org/content/10.1101/2023.07.14.549110v2). They're provided with the GPLv2 license, in the hope that others find them useful in their own projects, but also with the explicit intention of receiving improvements as others make them in their own contexts. Below you'll find a discussion of the general workflow from that paper, as well as some documentation of what each script does. Note that these command line tools all also have `argparse` I/O tools, so if you're wondering about how to use them in more detail calling them with no arguments, or with the `-h` or `--help` flags will print some help text that may be more rapidly useful than leafing through this document.

# Installation

These python scripts are designed to be used as command-line tools; some can also be imported into python interpreter sessions or notebooks, if one is interested in calling the functions from them independently. Thus, getting python environments set up with the correct modules present is the main installation task. Getting the parameterization program associated with the docking code you're interested in installed is the other major task.

## System preparation programs

The two preparation programs we've used are for autodock family programs, and generate PDBQT files. They are:

### OpenBabel

Get it from conda-forge:
`mamba install -c conda-forge openbabel`

Get it from apt (for debian-like systems):
`sudo apt install openbabel`

### `prepare_receptor` and `prepare_ligand`:

Because these programs _require_ python 2, and therefore don't play nicely with any other python tools written anytime recently, they are best installed from the tarball available from the [ADFR website](https://ccsb.scripps.edu/adfr/downloads/). You should follow the directions there to get them added to your path.

Note that you can view the [licenses](https://ccsb.scripps.edu/adfr/license/) these programs are being provided under on that same site. Because we only need `prepare_receptor` and `prepare_ligand`, these are available under LGPLv2 at time of writing.

## Docking programs

Right now our scripts are set up to either use [VINA](https://vina.scripps.edu/) or [SMINA](https://github.com/mwojcikowski/smina). In both cases the easiest way to get either set up is to install into a conda env for docking:
```
mamba create -n docking -c conda-forge smina vina
```

We have run into an issue with using VINA from `conda-forge` source; in particular, we found [there was a bug](https://github.com/ccsb-scripps/AutoDock-Vina/issues/90) that produced pretty significant variability in docking scores given the same pose, which is pretty bad for this application where the scores themselves are being used. We also noticed another bug whereby sometimes, seemingly at random across thousands of docking runs, a handful of poses will be given an extra 10-15 kcal/mol of favorability. There are several other energy function bugs on their issue tracker that this could be. Until the version on conda forge is greater than `1.2.4` we'd recommend compiling the `develop` branch of [VINA hosted on github](https://github.com/ccsb-scripps/AutoDock-Vina), following their install instructions for [building from source](https://autodock-vina.readthedocs.io/en/latest/installation.html#building-from-source), focusing on the python bindings since that's the only part of the install used by `docking_parallel.py`.

Note if you are trying to get the right compilers available you can install those from conda forge via the `compilers` repository. Here's an example install line:
```
mamba create -n built-vina -c conda-forge compilers boost-cpp swig numpy
```
Then you'd set up the python bindings:
```
conda activate built-vina
cd path/to/cloned/Autodock-Vina/build/python
rm -rf build dist *.egg-info
python setup.py build install
```

## Trajectory manipulation and postprocessing

To run the python tools associated with trajectory manipulation and post processing you'll need [`loos`](https://github.com/GrossfieldLab/loos), a trajectory handling library, and [`enspara`](https://github.com/bowman-lab/enspara), an MSM/utilities library. You'll also very likely need either [`pyemma`](http://www.emma-project.org/latest/), [`deeptime`](http://www.emma-project.org/latest/), or both.

A tidy way to install all of these into the same environment is to grab the libraries that are available on `conda-forge` and also the depeendencies that are needed for enspara, then clone enspara and build it using its `setup.py` script.

```
mamba create -n popshift -c conda-forge loos deeptime pyemma mdtraj cython mpi4py
conda activate popshift
git clone https://github.com/bowman-lab/enspara
cd enspara
python setup.py install
```

# Workflows

The scripts provided here are supposed to be modular enough that you can adjust them to your needs, or swap one part out in place of another part--we are certainly going to do this as our research on this framework evolves, so why shouldn't you? However, there are a few pathways that probably make general sense to understand, and the easiest way to begin to see what order things need to happen in is probably with some prescribed workflows. Find some below:

## Docking to an MSM

In the original PopShift MS, we considered how to apply the PopShift framework to the ensemble docking problem. The general workflow to actuate that study, with the tool(s) needed, follows:

1. Obtain a satisfactory MSM representing ligand-free simulations of your receptor.
2. Determine a selection that picks out your pocket of interest.
3. Pick and align frames from the MSM using `pick_align_frames.py` to ensure that all the residues that should be in your box are nicely co-aligned.
   - Check that things landed where you expected them to land using a visualizer like `pymol` and `draw_box.py`.
   - Noting that the alignment will recenter the coordinate system of all the models at the centroid of the atoms selected as the 'pocket' atoms, so boxes centered at or close to `0,0,0` are probably good choices in general.
4. Prepare receptors and ligands (using either OpenBabel or AutoDock Tools `prepare_receptor` and `prepare_ligand`), here you have some alternatives:
   - Use `prepare_ligand` and `prepare_receptor` on the command line directly with your favorite shell concurrency tool such as `xargs` or GNU `parallel`.
   - Alternatively, use `obabel` with the `-o pdbqt` flag.
5. Dock to each sample, saving the docking score in the file containing the ligand pose, using `docking_parallel.py`.
6. Extract scores from these output files and collate them into arrays using `extract_scores.py`
7. Use the score arrays to do PopShift calculations, with `popshift.py`
8. Plot stuff, compare to experiments, and generally *go nuts*, by reading the JSON that `popshift.py` outputs into your downstream scripts or notebooks.

Note that the docking here is done with either VINA or SMINA, so steps 4, 5, and 6 are somewhat particular to how those codes do that. You can of course color outside of the lines here, but you'll probably be adding some functionality to the scripts we have in order to do so.

### Frame picking using `pick_align_frames.py`

For the docking below to work as completely independent jobs, we extract frames representing each state, then save them into subdirectories where each subdirectory name corresponds to the MSM state index. Each file name consists of a trajectory index, then a frame index, separated by a dash. Here is an example of frames picked from a 25 state MSM built from five longish trajectories (obtained by calling `tree -P '*.pdb' -v receptor`) output abbreviated for clarity:

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
The commandline I used to create this was:

```
frame_picker=path/to/PopShift/pick_align_frames.py

python $frame_picker \
  --find-bonds 1.95,1.35\
  --make-receptor-sel-chain-A\
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

Here are the options and arguments, in order:
- `find-bonds`: find bonds using two different cutoffs, the first for heavy atoms and the second for hydrogen. Not necessary if bonds are present in the file used for 'model'. Better to have bonds written into PDBs so that the prep program doesn't need to build them.
- `make-receptor-sel-chain-A`: make the written PDBs have chain A identification. Some system prep tools require this.
- `assignments`: the assignments or discretized trajectory file, saved as a 'ragged array. Alternatively, an array of center indices saved as a `.npy`. To subsample randomly the full assignments file is obligatory, and **the order of the trajectories in the assignments file and the coordinate trajectory paths must match**.[^1]
- `align-resid-list`: this expects a text-file that `numpy.genfromtxt` can turn into a 1-D array of `int`s. These correspond to the `resid`s from the model that should be used in the iterative alignment. 
- `number-frames`: Select this many frames from each MSM bin.
- `mapping`: this corresponds to an enspara-style 'to' trim mapping; a dictionary containing the indexes of the original assignments trajectories as keys and the indexes they should be mapped to after ergodic trimming as values.
- `system_name`: The name you'd like to give to this particular extraction. Will become the directory name within which the `receptor` dir containing the tree of selected frames is written.
- `prot_masses.pdb`: this is the model file. For this particular command I was using a `prot_masses.pdb` generated by gromacs, which means it doesn't have connectivity which is why the example also contains `find-bonds`.
- `eq-probs.npy`: a numpy array containing equilibrium probabilities. Could also be a PyEMMA MSM saved as a pickle file, in which case `mapping` would not have been necessary. Used to communicate to the script how many states are active.
- `random`: the sampling mode requested here randomly picks frames from each msm bin without replacement. There are other options; consult the help statement. `centers` will switch the tool to operating on the provided coordinate trajectory as a centers trajectory.
- `!hydrogen`: the [LOOS selection string](https://grossfieldlab.github.io/loos/selections.html) to apply to the residues selected as the pocket for the purposes of alignment. This string means 'not hydrogen', or in other words 'all heavy atoms'.
- `traj_list.txt`: a text file containing paths to the coordinate trajectoies frames should be extracted from. A sequence of at least one paths may also be provided here, but as previously mentioned **the order of the trajectories in the assignments file and the coordinate trajectory paths must match**. Because of this I much prefer to create a list of trajectories and save it as a text file when I cluster, then pass that list around to subsequent tools. If operating in centers mode then this argument should be a path to the trajectory of centers.

Calling `pick_align_frames.py --help` will provide more information about alternative strategies for picking frames (for example, using just cluster centers, or picking particular numbers of frames from each bin). It's also worth noting that pick_align_frames can be imported into other python sessions, so if you want to rip random frames into sub-trajectories for other reasons this is possible. Finally, note that subsample trajectories that are not PDBs can also be written to disk if desired. This can be nice if trying to use normal trajectory analysis tools on the selected frames, because reading in piles of PDBs instead of one DCD per MSM bin can be annoying and slow.

[^1]: The order of `*.xtc`, python's `glob.glob('*.xtc')` python's `pathlib.Path.glob('*.xtc')` and `find . -name '*.xtc'` will not necessarily be consistent across runs. The first of these results in ASCII-betical ordering (meaning that `traj-11.xtc` will come _before_ `traj-2.xtc` for bash and derivative shells). The other three, if unsorted, will come in the file-system's order. Importantly, **the filesystem order can change** if the contents of the filesystem change, which is why many references point out that the returns of such commands--when unsorted explicitly--are arbitrary. This is an easy way for one to shoot oneself in the foot.
# Docking scripts

## On preparing receptors and ligands for docking

We prepare each of the receptor conformations as a completely independent system, so each is represented by its own PDBQT file and can thus be run as an isolated job on whatever hardware is available. This leaves some things to be desired, but as a first pass it's worked for us.

In general these scripts are set up to play nicely with one another by assuming that the directory structure of the intermediate files will be maintained--this isn't a hard requirement, but it does keep things convenient. Because we are using other command-line tools to prepare ligands and receptors for docking, we are really just expecting that the prepared files will be in the same place as the files they were prepped from. There are a lot of ways to call a command-line tool on a file such that it writes a new file with the same path but a different file extension. Here are several I have used in the past:

Using GNU parallel:

```
# use find to get the paths rel. to CWD;
# write to file to prevent arg-list too long
find frame_picking_system_dir/receptor/ -name '*.pdb' > receptor_paths.txt
# Don't use too many jobs; I/O intensive task will be bound by disk perf.
parallel -j 8 obabel -opdbqt {} -O {}qt :::: receptor_paths.txt
```

Or with `prepare_receptor`:

```
parallel -j 8 prepare_receptor -r {} -o {}qt :::: receptor_paths.txt
```

What's going on here? Scrutinizing the docs for [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) suggests that parallel will be running the commandline on each of the elements in the file when in `::::` input mode. It substitutes each line into the 'command line' you provide as an argument wherever the curly braces are, then runs that as its own subprocess.

Using [xargs](https://manpages.org/xargs) as an alternative to parallel:

```
find frame_picking_system_dir/receptor/ -name '*.pdb' | xargs -I% -n 1 -P 8 obabel -opdbqt % -O %qt
```

If only dealing with a few frame-picking runs (or maybe just one) for a relatively parsimonious MSM (one with few states), we can afford to do this sequentially (nicer if you're getting errors, as well). In this case, by far the easiest thing to do is write a simple loop in the shell of your choosing. Here's a while loop that takes advantage of the redirect-to-file from `find` given in the `parallel` example:

```
while read receptor_path; do
    echo $receptor_path
    obabel -opdbqt $receptor_path -O ${receptor_path}qt
done < receptor_paths.txt
```

### Preparing Ligands

Ligands can also be prepped by command-line tools; the docking script `docking_parallel.py` expects you to provide the ligands as either a text-file of ligand paths or as a series of command-line arguments. As such, you can organize these files how you see fit. I usually put them in their own directory, then prep 'everything' with a certain file extension within that directory. So for example:

```
# Assume there's a directory called ligands with .sdf or .mol2 files in it
parallel -j 8 obabel {} -h -isdf -opdbqt -O{.}.pdbqt ::: ligands/*.sdf
```
NOTE: the `-h` flag adds hydrogens. Don't throw it if you don't want that.

Alternatively, if you have the ADFR suite installed and would like to use `prepare_ligand` instead of openbabel (note you need `.mol2` or `.pdb` files for `prepare_ligand`, as of this writing):

```
parallel -j 8 prepare_ligand -l {} -Ahydrogens -o{.}.pdbqt ::: ligands/*.mol2
```

## Docking using Jug

Again here there are many ways to organize docking so long as the docking code is available on the command-line. We like to use Jug, because it integrates well with our workload manager (slurm). We have also used it with lsf. Following is an example sbatch script, with comments about the various options for docking paprallel.
```
#!/bin/bash
#SBATCH --array=0-99
#SBATCH -n 1
#SBATCH -p all
#SBATCH -o dockouts/%A-%a.out
#SBATCH -e dockouts/%A-%a.out

receptor_dir=system_name_from_framepicking
jug execute -- path/to/PopShift/docking_parallel.py\
 -d smina\
 -e 32\
 $receptor_dir/receptor\
 $receptor_dir/'12xsmina'\
 '0,0,0'\
 '12,12,12'\
 all-ligands/*.pdbqt

```
The SBATCH pragma lines are just asking for a 100 element array job where each job is a one-core job with default everything on the 'all' partition. Obviously you can adjust that if you'd like, but the tool expects to just run single-threaded docking jobs, so there's no point to increasing the per-job core count with the script as written. The command lines are as follows:
- `jug execute --`: saying this invokes jug instead of the python interpreter, issuing the execute command. The `--` tells the jug program that the part of the command line it should read is concluded; if this is omitted and the command line contains long-option flags that should be going to `docking_parallel.py` jug will interpret them as being for it instead.
- `d` is the 'docking-algorithm' flag, and picks from the available docking options--at present either `smina` or `vina`.
- `e` is the 'exhaustiveness' flag; for vina, the exhaustiveness flag is the number of MC runs to perform.
- `$receptor_dir/receptor`: the directory containing the receptor conformations to dock to.
- `$receptor_dir/12xsmina`: the nickname to give the docking run. It's preferable to keep all the docking runs organized in the same directory to their `receptor` dir for ease of later scripting. 
- `0,0,0`: the coordinates of the center, separated by commas. Note that if you used `pick_align_frames.py` to get receptor structures to dock to, 0,0,0 is a relatively good choice for a center because the coordinates get zeroed by the frame 
- `12,12,12`: the edge lengths of the x, y, and z edges of the docking box.
- `all-ligands/*.pdbqt`: the last N arguments should be prepped ligand PDBQTs. I bunged the prepped PDBQTs of the ligands into a directory called `all_ligands`, then provided them to the command-line as a glob. Note that which order the ligands get docked in is not in general important.

### Output

The output from docking will be written into a directory structure that is within a directory given to the nickname for the docking run, with sub-directories using the stem of the file-names of the ligand PDBQTs, and the directory structure beneath that containing numbered directories for each MSM bin, and each docked pose containing _only_ the ligand coordinates will have a file-name matching the receptor structure to which it was docked. So for the above frame-picking example, supposing the ligand files docked in were `benzene.pdbqt`, `toluene.pdbqt`, and `p-xylene`:

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

 ### Issuing Jug commands
There are other [jug commands](https://jug.readthedocs.io/en/latest/subcommands.html). In order to issue them for a jug file that has command line arguments, one must make sure that the same command-line arguments are provided. One must also be using the same versions of python and jug, so its important to have the environment the dock job was launched from active. Thus, often the easiest thing to do is make the sbatch into a script that takes an argument, and inserts that argument in place of `execute`:

```
#!/bin/bash

receptor_dir=system_name_from_framepicking
jug $1 -- path/to/PopShift/docking_parallel.py\
 -d smina\
 -e 32\
 $receptor_dir/receptor\
 $receptor_dir/'12xsmina'\
 '0,0,0'\
 '12,12,12'\
 all-ligands/*.pdbqt
```
Call this script `issue-jug-cmds-dp.sh` or similar. I could then use this script to obtain jug status reports, or to invalidate the results of a calculation so that I could rerun it:
```
issue-jug-cmds-dp.sh status
issue-jug-cmds-dp.sh 'invalidate dock_smina'
issue-jug-cmds-dp.sh cleanup
```

## Analyze results with PopShift framework

Now we are ready to do something with the per-state affinity and pose estimates we've made. We'll extract the scores, then use them to make affinity estimates.

### Extract scores

Extract scores with `extract_scores.py`. This will pull the values of docking scores out of the various pose files and save them to an [`enspara`](https://enspara.readthedocs.io/en/latest/installation.html) ragged array in h5 format, since it is not in general true that the number of frames extracted per bin is always the same. This is a pretty simple tool; here is an example command line:

```
python path/to/PopShift/extract_scores.py -n 8 -t smina 12xsmina
```
The commandline works as follows:
- `n`: the number of cores to use; this job is pretty I/O bound, but I have found that using a handful of cores can speed it up. 
- `t`: the type of output file the script can read. Will match the `d` option from `docking_parallel.py`.
- `12xsmina`: the docking run nickname we would like to run the extraction out of.
  
By default, the resulting score extracts will be written to `extracted_scores`, in the same directory as the docking run is in. It takes the path to that directory, then goes up one level. Within `extracted_scores` there will be a dir for each docking run you've extracted, and it will be named to match the directory name for the docking run. There will be an .h5 file within this directory for each ligand name. The elements will correspond to the scores extracted from the pose files.

### Calculate macroscopic binding constants and reweighted state probabilities

To obtain binding free energy estimates and reweighted state populations, use `popshift.py`. This program was designed to either be a commandline tool that will perform these analyses on files arranged in the expected way, or as a small module with functions to perform the desired operations, if working within a notebook is desired. Here's an example of command-line usage, but as with all the other tools this one has broader summary documentation provided by the output of the `--help` flag.

```
extracts=extracted_scores/12xsmina
python path/to/PopShift/popshift.py -n 8\
  --out $extracts/binding-calx\
  --reweighted-eq 1\
  bin-samples\
  path/to/eq-probs.npy\
  $extracts/*.h5
```


## Dependencies
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
    - > conda install -c conda-forge smina
4. popshift.py
   - Install enspara into your python environment, so that `popshift.py` can import it when run
   - If using PyEMMA discretized trajectories, then it must be installed:
     - > conda install -c conda-forge pyemma
   - Usage can be had by saying `python popshift.py` with no arguments. Help can be had with the `--help` flag.
5. LOOS (for `pick_align_frames.py`, `draw_box.py`, `add_bonds_posthoc.py`, and `rmsd_receptor_ligand.py`)
    - > conda install -c conda-forge loos
    - Can also be compiled following the instructions on the [loos github](https://github.com/GrossfieldLab/loos)
## Conda command for dependencies
### minimal
   > mamba install -c conda-forge -c insilichem autodocktools-prepare smina loos cython mdtraj 

- these last two will be for installing enspara, which needs to be done by cloning the repository and using the `setup.py` script with the correct conda environment activated.

- Note: you still need to do compilation steps for enspara and (if desired) autodock VINA.   

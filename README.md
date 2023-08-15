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
   - `prep_parallel.py` uses python's `multiprocessing` module to call `prepare_receptor`. `prepare_ligand`, because it is faster, was not parallelized in this way.
   - Alternatively, use `prepare_ligand` and `prepare_receptor` on the command line directly with your favorite shell concurrency tool such as `xargs` or GNU `parallel`.
   - Alternatively, use `obabel` with the `-o pdbqt` flag.
5. Dock to each sample, saving the docking score in the file containing the ligand pose, using `docking_parallel.py`.
6. Extract scores from these output files and collate them into arrays using `extract_scores.py`
7. Use the score arrays to do PopShift calculations, with `popshift.py`
8. Plot stuff, compare to experiments, and generally *go nuts*, by reading the JSON that `popshift.py` outputs into your downstream scripts or notebooks.

Note that the docking here is done with either VINA or SMINA, so steps 4, 5, and 6 are somewhat particular to how those codes do that. You can of course color outside of the lines here, but you'll probably be adding some functionality to the scripts we have in order to do so.

### Frame picking using `pick_align_frames.py`

For the docking below to work as completely independent jobs, we extract frames representing each state, then save them into subdirectories where each subdirectory name corresponds to the MSM state index. Each file name consists of a trajectory index, then a frame index, separated by a dash. Here is an example of frames picked from a 25 state MSM built from five longish trajectories (obtained by calling `tree -v receptor`) output abbreviated for clarity:

```
receptor/
├── 0
│   ├── 0-55257.pdb
│   ├── 0-55257.pdbqt
│   ├── 0-56428.pdb
│   ├── 0-56428.pdbqt
│   ├── 1-43334.pdb
│   └── 1-43334.pdbqt
├── 1
│   ├── 2-105110.pdb
│   ├── 2-105110.pdbqt
│   ├── 2-175513.pdb
│   ├── 2-175513.pdbqt
│   ├── 3-48068.pdb
│   └── 3-48068.pdbqt
├── 2
│   ├── 0-157652.pdb
│   ├── 0-157652.pdbqt
│   ├── 0-158908.pdb
│   ├── 0-158908.pdbqt
│   ├── 0-191040.pdb
│   └── 0-191040.pdbqt
...
└── 24
    ├── 4-24060.pdb
    ├── 4-24060.pdbqt
    ├── 4-186351.pdb
    ├── 4-186351.pdbqt
    ├── 4-189200.pdb
    └── 4-189200.pdbqt

25 directories, 150 files
```

### Preparing Receptors

# Docking scripts

## On preparing receptors and ligands for docking

We prepare each of the receptor conformations as a completely independent system, so it gets its own PDBQT file and can thus be run as an isolated job on whatever hardware is available. This leaves some things to be desired, but as a first pass it's worked for us.

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
parallel -j 8 prepare_receptor -r {} -o {}.pdbqt :::: receptor_paths.txt
```

What's going on here? Scrutinizing the docs for [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) suggests that parallel will be running the commandline on each of the elements in the file when in `::::` input mode. It substitutes each line into the 'command line' you provide as an argument wherever the curly braces are, then runs that as its own subprocess.

Using xargs as an alternative to parallel:

```
find frame_picking_system_dir/receptor/ -name '*.pdb' | xargs -I% -n 1 -P 8 obabel -opdbqt % -O %qt
```

If only dealing with a few frame-picking runs (or maybe just one) for a relatively parsimonious MSM (one with few states), we can afford to do this sequentially (nicer if you're getting errors, as well). In this case, by far the easiest thing to do is write a simple loop in the shell of your chosing. Here's a while loop that takes advantage of the redirect-to-file from `find` given in the `parallel` example:

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

Again here there are many ways to organize docking so long as the docking code is available on the command-line. We like to use Jug, because it integrates well with our workload manager (slurm). We have also used it with lsf.


## Analyze results with MSM docking framework
### Extract scores
- Extract scores with `extract_scores.py`
- help available with `--help`
- Should be provided each set of docked poses for each ligand as a group.
### Calculate macroscopic binding constants and reweighted state probabilities
- `popshift.py` for both binding free energies and reweighted populations.


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

# Introduction

This is a repository containing some command-line tools for implementing the PopShift framework for docking. The purpose behind these tools are to allow users both to make alterations to these workflows for their own use, and also to apply them in similar fashion to the use they were put toward in the original [PopShift manuscript](https://www.biorxiv.org/content/10.1101/2023.07.14.549110v2). They're provided with the GPLv2 license, in the hope that others find them useful in their own projects, but also with the explicit intention of receiving improvements as others make them in their own contexts. Below you'll find a discussion of the general workflow from that paper, as well as some documentation of what each script does. Note that these command line tools all also have `argparse` I/O tools, so if you're wondering about how to use them in more detail calling them with no arguments, or with the `-h` or `--help` flags will print some help text that may be more rapidly useful than leafing through this document.

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

### On preparing ligands

In general these scripts are set up to play nicely with one another by assuming that the directory structure of the intermediate files will be maintained--this isn't a hard requirement, but it does keep things convenient. Because we are using other command-line tools to prepare ligands and receptors for docking, we are really just expecting that the prepared files will be in the same place as the files they were prepped from. There are a lot of ways to call a command-line tool on a file such that it writes a new file with the same path but a different file extension. Here are several I like/have used in the past:

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

What's going on here? Scrutinizing the docs for [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) suggests that parallel will be running the commandline 

Sometimes we don't need this done at breakneck, and can afford to simply do this sequentially. In this case, by far the easiest thing to do is write a simple loop in the shell of your chosing. Here's a while loop that takes advantage of the file-redirect from find we used earlier:
```
while read receptor_path; do
    obabel -opdbqt $receptor_path -O ${receptor_path}qt
done < receptor_paths.txt
``````

# Docking scripts

## Prepare receptor and ligands
prep_parallel.py --> use to prepare receptor files and ligands

    python prep_parallel.py --help # for options on how to run

prepare_ligand.py --> prep_parallel does this now, so this is not needed  

## Dock
docking_parallel.py     	--> use for docking in parallel on a single node using vina or smina scoring function  

    python docking_parallel.py --help # for options on how to run

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
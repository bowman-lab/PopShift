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
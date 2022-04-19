# Docking scripts

## Prepare receptor and ligands
prep_parallel.py --> use to prepare receptor files and ligands

    python prep_parallel.py --help # for options on how to run

prepare_ligand.py --> prep_parallel does this now, so this is not needed  

## Dock
docking_parallel.py     	--> use for docking in parallel on a single node using vina or smina scoring function  

    python docking_parallel.py --help # for options on how to run

## Analyze results with MSM docking framework
### extract scores
extract-scores.sh 

## Dependencies
1. ADFR suite
    - > https://ccsb.scripps.edu/adfr/downloads/  
    cp prepare_ligand and prepare_receptor into $PATH
2. Autodock Vina python bindings
    - > git clone https://github.com/ccsb-scripps/AutoDock-Vina  
    cd AutoDock-Vina/build/python  
    conda install -c conda-forge numpy boost-cpp swig  
    rm -rf build dist *.egg-info (to clean previous installation)  
    python setup.py build install  
3. Smina
    - > conda install -c conda-forge smina
4. Antechamber from AmberTools 
    - > conda install -c conda-forge ambertools
5. extract-scores.sh
    - parallelized using GNU parallel (say `which parallel` at the command line to see if you have it). 
    - If you don't have it, you can install using `conda install -c conda-forge parallel`.
    - If using **bash**, the following one-liner will check and install for you (assuming you're in the right env):
    - > type parallel 2>/dev/null || { echo >&2 "Needed to install parallel using conda-forge."; conda install -c conda-forge parallel; }
6. msm_binder.py
   - Install enspara into your python environment, so that `msm_binder.py` can import it when run
   - If using PyEMMA discretized trajectories, then it must be installed:
     - > conda install -c conda-forge pyemma
   - Usage can be had by saying `python msm_binder.py` with no arguments. Help can be had with the `--help` flag.
### Install conda dependencies in one line (better results; conda less likely to get confused):
   >mamba install -c conda-forge numpy boost-cpp swig smina ambertools parallel pyemma
- Note: you still need to do compilation steps for enspara and autodock vina, most likely. You can do this all as one line first, though.   
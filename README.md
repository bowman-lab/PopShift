# Docking scripts

## Prepare receptor and ligands
prep_parallel.py --> use to prepare receptor files and ligands

    python prep_parallel.py --help # for options on how to run

prepare_ligand.py --> prep_parallel does this now, so this is not needed  

## Dock
docking_parallel.py     	--> use for docking in parallel on a single node using vina or smina scoring function  

    python docking_parallel.py --help # for options on how to run

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
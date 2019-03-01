#!/bin/bash   
source /ChemistryData/jgwdk/lib/Anaconda/bin/activate my-rdkit-env
source /ChemistryData/jgwdk/lib/Anaconda/bin/activate my-molvs-env
python cleanPredSmiles.py $1

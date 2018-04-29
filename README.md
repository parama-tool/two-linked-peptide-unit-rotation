# two-linked-peptide-unit-rotation

This repository contains the python script which can be used to "recreate" the Ramachandran map from scratch
It also contains the dataset of coordinates of 1789 two-linked peptides units (peptide_unit_coordinates.txt), identified from small molecule crystal structures in Cambridge Structural Database. 

To run the script type the following in command line:
  python peptide_rotation.py <file_name>

******************************************************************************************************
PLEASE NOTE:
1. Make sure the input file is in PDB format and present in the same directory as the python script. Read through the comments in the code for more information on input file requirements
2. Two output files will be created - the shrt_conts.txt file contains the short contact information for every partially allowed or disallowed phi-psi combination. The rmap.csv file contains information on whether a phi-psi combination is fully allowed, partially allowed or disallowed

*******************************************************************************************************

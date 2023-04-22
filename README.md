# hps2gmx_topol_converter
[![GitHub](https://img.shields.io/github/license/sscake/hps2gmx_topol_converter)](LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/sscake/hps2gmx_topol_converter)](https://www.python.org/)  
This program aims to generate the topology file in GROMACS format to use GROMACS built-in analysis program after the implementation of HydroPathy Scale (HPS) simulation.  

The HPS parameters of amino acids used in this program are from [here](https://doi.org/10.1371/journal.pcbi.1005941), and the ones of nucleotides are from [here](https://doi.org/10.1093/nar/gkaa1099). 

# Preparation
This program have to read the HPS topololgy information to generate the GROMACS topology file, thus the HPS topology information should be provided in the following format.  

```
## Starting of the topology information ##  
## Any information about the simulation. ##  
types   <type_number>    # type_number is the total number of the bead types, e.g., 23  
<a list of all the types>    
# Amino acids are represented by uppercase one letter. If the amino acid is at 
# N/C-terminal, "N/C" have to be appended after the one-letter representation. In 
# addition, nucleotides of RNA are represented by lowercase one letter, and "5" 
# should be appended after the nucleotide at 5'-terminal. For instance, 
# Y  MC  NN  S  R  P  W  Q  G  D  F  E  A  K  I  N  L  M  g  u  c  a  g5  

atom   <atom_number>  
# index name mass  charge  
<list of all atoms (beads), zero-based indexing>  
# e.g., 0   NN 114.10    1  
#       1    R 156.20    1  
.  
.  
.  
  
pbonds   <total bond number of amino acids>  
rbonds   <total bond number of nucleotides>  
bonds    <total bond number>  
<list of all bonds between two atoms (beads), zero-based indexing>  
# e.g., 0       1  
#       1       2  
.  
.  
.  
  
## End of the topology information ##  
```
  
  
# Dependency
This program has dependency on "cython" package, thus it should be installed before using this program.   
*pip install cython* can be used to install it. If you have installed conda, you can also use the command: *conda install cython* .   

# Installation
Turn to your working directory, i.e., *cd ...*  
Give the *setup* file executable permission, i.e., *chmod u+x setup* .  
After that, it can be easily installed by typing in your console: *./setup build_ext --inplace* .  

# Usage
Give the *run_converter* file executable permission as well, i.e., *chmod u+x run_converter* .  
The uasge is also very convenient. Just type *./run_converter -hps /input/path/to/your/hps_topology/file -topol /output/path/to/your/gromacs/topol/file* .  
  
If your topology information filename is "hps_topol.txt" and the GROMACS topology filename is expected to be "topol.top", just run *./run_converter* !  
  
Hope you enjoy it!   

If you encounter any problem when using it, feel free to contact me!  
E-mail: xsliu16@fudan.edu.cn

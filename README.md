# shaped_penetrant

This is the GitHub repository for studying shaped penetrant diffusion in highly crosslinked networks. For more detailed information, please refer to J. Chem. Phys. 160, 114905 (2024) https://doi.org/10.1063/5.0197140

# Create initial structures

## 1. Compile with 

```
g++ -g -o in-cfg-nonsph-cpp.out init_config_nw_nonsph.cpp
```

## 2. Create structures

### use the **mk-in-cfg-rod.sh** direclty under the same directory

```
bash mk-in-cfg-rod.sh Nrod AR len_pene Dp
```

### create one initial configuration at a time

Note: Make sure the folder named **f_${fcross[i]}/initcfg-Nrod${Nrod}-AR${AR}-Dp${Dp}-lp${len_pene}** has been created before running the following 

Run the code from the same directory using the executable with the following parameters: 

    N_chain: number of chains
    N_mono: umber of beads in each chain  
    N_cross: number of reactive bead on each chain
    N_dangle: number of dangling bead on each reative bead
    Nrod: number of rod-like penetrants
    AR: how many bead a penetrants is composed of
    len_pene: chain length between penetrant beads
    Dp: penetrant bead diameter

```
./in-cfg-nonsph-cpp.out 60 30 6 0 3 2 1.587 1.587
```

# Run simulations on LAMMPS

## Create crosslinked networks at T=1

## Study diffusion at target T

# Raw data for figures in the main article and supporting information

# Authors
Tsai-Wei Lin and Charles E. Sing

# Funding Acknowledgements

This research was supported by the U.S. Department of Energy, Office of Basic Energy Sciences, Division of Materials Sciences and Engineering (Award No. DE-SC0020858), through the Materials Research Laboratory at the University of Illinois at Urbana-Champaign. Helpful discussions with Kenneth Schweizer and Baicheng Mei are gratefully acknowledged.


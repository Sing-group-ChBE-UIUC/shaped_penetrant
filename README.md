# shaped_penetrant

This is the GitHub repository for studying shaped penetrant diffusion in highly crosslinked networks. For more detailed information, please refer to J. Chem. Phys. 160, 114905 (2024) https://doi.org/10.1063/5.0197140

# Create initial configuration

## 1. Compile with 

```
g++ -g -o in-cfg-nonsph-cpp.out init_config_nw_nonsph.cpp
```

## 2. Run the code from the same directory using the executable with the following parameters: 
    N_chain = atoi(argv[1]);    // number of chains
    N_mono = atoi(argv[2]);     // number of beads in each chain  
    N_cross = atoi(argv[3]);    // number of reactive bead on each chain
    N_dangle = atoi(argv[4]);   // number of dangling bead on each reative bead
    Nrod = atoi(argv[5]);       // number of rod-like penetrants
    AR = atoi(argv[6]);         // how many bead a penetrants is composed of
    len_pene = atof(argv[7]);   // chain length between penetrant beads
    Dp = atof(argv[8]);         // penetrant bead diameter
    
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


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

Use Dp=1.587, Nrod=3, Enp=1.0, fcross = 0.333 as an example

## 1. Copy the initial configuration to each folder

## 2. Create crosslinked networks at T=1

Under the folder **Dp1.587-Np3-Enp1.0/f_0.333/T1.0-387279973-xln**, run
```
lmp_serial -in p-soft-T1
```
and then run 

```
lmp_serial -in p-lj-T1-xln
```

## 3. Study diffusion at target T
Under the folder **Dp1.587-Np3-Enp1.0/f_0.333/T0.7-387279973-nvt**, run in the following order

#### i. copy the **re-nw-T1.0** to folder **T{Ttarget}-{seednumber}-nvt**

#### ii. cool the system to target temperature

```
lmp_serial -in p-lj-T0.7-cool
```
#### iii. short npt run to get averaged V

```
lmp_serial -in p-lj-T0.7-short-npt
```

#### iv. production run at NVT

```
lmp_serial -in p-lj-T0.7-nvt
```

#### v. short NVT run to get short time data (optional)

```
lmp_serial -in p-lj-T0.7-isf
```

# Raw data for figures in the main article and supporting information

# Authors
Tsai-Wei Lin and Charles E. Sing

# Funding Acknowledgements

This research was supported by the U.S. Department of Energy, Office of Basic Energy Sciences, Division of Materials Sciences and Engineering (Award No. DE-SC0020858), through the Materials Research Laboratory at the University of Illinois at Urbana-Champaign. Helpful discussions with Kenneth Schweizer and Baicheng Mei are gratefully acknowledged.


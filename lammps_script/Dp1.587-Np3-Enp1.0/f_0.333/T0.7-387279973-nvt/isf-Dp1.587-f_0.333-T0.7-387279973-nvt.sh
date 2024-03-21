#!/usr/bin/env bash

#Shell
#!/bin/bash

#Combine STDERR and STDOUT
#$ -j y

#Send mail on abort(a), begin(b), and end(e)
#$ -M twlin2@illinois.edu

#Select queue
#$ -q xeon1.q,xeon2.q,serial.q,infiniband.q

#Run job from the directory from which it was submitted
#$ -cwd

#Set priority
#$ -p -10

#Job name
#$ -N isf-Dp1.587-f_0.333-T0.7-387279973-nvt

#Run the executable
module load gcc/9.3.0
~/applications/lammps-stable_3Mar2020/src/lmp_serial -in p-lj-T0.7-isf

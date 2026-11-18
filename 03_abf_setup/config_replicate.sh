#!/bin/bash
# Generates one NAMD config per ABF window. The config layout is adapted from
# the NAMD tutorial configuration files (Theoretical and Computational
# Biophysics Group, University of Illinois Urbana-Champaign).

base_dir="$(pwd)"   # run from the system directory

for ((i=2; i<=2; i++))
do
    filename="win${i}.conf"
    inputname="Can_gold-${i}"
    restartname="Can_gold-${i}"

    file_path="${base_dir}/${filename}"

    cat > "${file_path}" <<EOF
#############################################################
## JOB DESCRIPTION
#############################################################

# ABF of Normal methane gold NP, window ${i}

#############################################################
## ADJUSTABLE PARAMETERS
#############################################################

structure          ../gold_can_membrane.psf
coordinates        ../gold_can_membrane.pdb

exec mkdir -p output

set temperature    310
set outputname     output/win${i}
set inputname      ${inputname}
set restartname    ${restartname}

if {1} {
binCoordinates     \$restartname.restart.coor
extendedSystem     \$restartname.restart.xsc
}

firsttimestep 0

#############################################################
## SIMULATION PARAMETERS
#############################################################

paraTypeCharmm	    on

parameters              /path/to/toppar/par_all36_prot.prm
parameters              /path/to/toppar/toppar_all36_carb_glycolipid.str
parameters              /path/to/toppar/toppar_all36_lipid_sphingo.str
parameters              /path/to/toppar/toppar_all36_lipid_detergent.str
parameters              /path/to/toppar/toppar_all36_lipid_prot.str
parameters              /path/to/toppar/par_all36_na.prm
parameters              /path/to/toppar/par_all36_lipid.prm
parameters              /path/to/toppar/par_all36_carb.prm
parameters              /path/to/toppar/par_all36_cgenff.prm
parameters              /path/to/toppar/par_interface.prm   ;#when added gold
parameters              /path/to/toppar/toppar_all36_lipid_cholesterol.str  ;#when added cholesterol
parameters 		/path/to/toppar/toppar_water_ions_cufix.str
parameters 		/path/to/toppar/toppar_all36_dphpc.str 


temperature         \$temperature

rigidBonds          all
exclude             scaled1-4
1-4scaling          1.0
switching           on
switchdist          8
cutoff              10
pairlistdist        12
margin               5

timestep            2.0
nonbondedFreq       1
fullElectFrequency  3
stepspercycle       12

langevin            on
langevinDamping     1
langevinTemp        \$temperature
langevinHydrogen    off

wrapAll             off

PME                 yes
PMEGridSpacing      1.0

langevinPiston        on
langevinPistonTarget  1.01325
langevinPistonPeriod  200.
langevinPistonDecay   100.
langevinPistonTemp    \$temperature

useGroupPressure      yes
useFlexibleCell       yes
useConstantRatio      yes

#############################################################
## OUTPUT
#############################################################

outputName          \$outputname

restartfreq         9600
dcdfreq             9600000
xstFreq             9600
outputEnergies      9600
outputPressure      9600

#############################################################
## EXTRA PARAMETERS
#############################################################

colvars              on
colvarsConfig        \${inputname}.in

#############################################################
## EXECUTION SCRIPT
#############################################################

minimize 48000

reinitvels \$temperature

run 20000004 

EOF

done

echo "Files created successfully."

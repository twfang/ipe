# csh shell commands to define the environment variables and aliases to run IPE 
module purge
module load intel
module load mvapich2
setenv COMPILER ifort
echo $COMPILER
setenv machine jet_$COMPILER
echo $machine

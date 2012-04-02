# Bash shell commands to define the environment variables and aliases to run IPE 
module load intel
export COMPILER=ifort
echo $COMPILER
export machine=zeus_$COMPILER
echo $machine

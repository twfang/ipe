# Bash shell commands to define the environment variables and aliases to run IPE 
#export COMPILER=pgf95
#export COMPILER=lf648    # lf64 v80
export COMPILER=ifort
echo $COMPILER
export machine=jet_$COMPILER
echo $machine

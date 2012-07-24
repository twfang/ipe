# Bash shell commands to define the environment variables and aliases to run IPE 
module switch intel lahey
export COMPILER=lf95
echo $COMPILER
export machine=jet_$COMPILER
echo $machine

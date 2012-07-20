# csh shell commands to define the environment variables and aliases to run IPE 
module switch intel lahey
setenv COMPILER lf95
echo $COMPILER
setenv machine jet_$COMPILER
echo $machine

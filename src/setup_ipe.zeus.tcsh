# Bash shell commands to define the environment variables and aliases to run IPE 
module load intel
setenv COMPILER ifort
echo $COMPILER
setenv machine zeus_$COMPILER
echo $machine

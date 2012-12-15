# csh shell commands to define the environment variables and aliases to run IPE 
module purge
module load lahey
#module load mpt
setenv COMPILER lf95
echo $COMPILER
setenv machine zeus_$COMPILER
echo $machine

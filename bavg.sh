#/bin/bash

MOAB_PATH=/root/opt/moab/
MOAB_LIBRARY=$MOAB_PATH"lib"
MOAB_INCLUDE=$MOAB_PATH"include"

DAGMC_PATH=/root/opt/dagmc/
DAGMC_LIBRARY=$DAGMC_PATH"lib"
DAGMC_INCLUDE=$DAGMC_PATH"include"
echo $MOAB_LIBRARY

g++  -std=c++11 avg.cpp -g -I$MOAB_INCLUDE -L$MOAB_LIBRARY -lMOAB -I$DAGMC_INCLUDE -L$DAGMC_LIBRARY -ldagmc -lboost_program_options -o avg


export DPSIM_ROOT=$(cd $(dirname ${BASH_SOURCE[0]}); pwd)

export PATH=$DPSIM_ROOT/inc:${PATH} 
export LD_LIBRARY_PATH=$DPSIM_ROOT/build/lib:${LD_LIBRARY_PATH} 


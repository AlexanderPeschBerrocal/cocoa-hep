export CURRENTDIR=$(pwd)

export LCG_109_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
source ${LCG_109_VIEW}/setup.sh

# Keep the legacy variables used by COCOA's CMake files while letting the
# LCG view provide the full, internally consistent runtime environment.
export GCC_HOME=/cvmfs/sft.cern.ch/lcg/releases/LCG_109/gcc/14.3.0/x86_64-el9
export PATH=${GCC_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${GCC_HOME}/lib64:${GCC_HOME}/lib:${LD_LIBRARY_PATH}
export CC=${GCC_HOME}/bin/gcc
export CXX=${GCC_HOME}/bin/g++
export ROOT_DIR=${ROOTSYS}/cmake
export GEANT4_HOME=${LCG_109_VIEW}
export GEANT4_DIR=${GEANT4_HOME}/lib64/cmake/Geant4
export CLHEP_HOME=${LCG_109_VIEW}
export CLHEP_DIR=${CLHEP_HOME}
export VECGEOM_HOME=${LCG_109_VIEW}
export VECGEOM_DIR=${VECGEOM_HOME}/lib64/cmake/VecGeom
export XERCESC_HOME=${LCG_109_VIEW}
export XERCESC_DIR=${XERCESC_HOME}
export VC_HOME=${LCG_109_VIEW}
export VC_DIR=${VC_HOME}/lib/cmake/Vc
export VDT_HOME=${LCG_109_VIEW}
export VDT_INCLUDE_DIR=${VDT_HOME}/include
export VDT_LIBRARY=${VDT_HOME}/lib/libvdt.so
export TBB_HOME=${LCG_109_VIEW}
export TBB_DIR=${TBB_HOME}/lib64/cmake/TBB
export HEPMC_HOME=${LCG_109_VIEW}
export JSONCPP_HOME=${LCG_109_VIEW}
export PYTHIA8_HOME=${PYTHIA8:-${LCG_109_VIEW}}
export FASTJET_HOME=$(fastjet-config --prefix)
export Qt5_DIR=${LCG_109_VIEW}/lib/cmake/Qt5

eval "$(geant4-config --sh)"
export PYTHIA8DATA=${PYTHIA8_HOME}/share/Pythia8/xmldoc

cd $CURRENTDIR

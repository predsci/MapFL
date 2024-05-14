#!/bin/bash
###################################################################
# Please set the location of HDF5 include/library files and
#  the linker flags to match your installed version.
#
# Note! The HDF5 library needs to have been compiled with
#  the same compiler being used here and is loaded in the run-time
#  environment (e.g. LD_LIBRARY_PATH).
#################################################################

# Location of local hdf5 installed with same compiler being used for POT3D:
HDF5_INCLUDE_DIR="/usr/local/include/"
HDF5_LIB_DIR="/usr/local/lib/"
# Fortran HDF5 library flags (these can be version dependent):
HDF5_LIB_FLAGS="-lhdf5_fortran -lhdf5_hl_fortran -lhdf5 -lhdf5_hl"

###########################################################################
# Please set the compile and flags based on your compiler and hardware setup.
###########################################################################

FC="gfortran"
FFLAGS="-O3 -march=native -fopenmp"

###########################################################################
###########################################################################
###########################################################################

MAPFL_HOME=$PWD

cd ${MAPFL_HOME}/src
echo "Making copy of Makefile..."
cp Makefile.template Makefile
echo "Modifying Makefile to chosen flags..."
sed -i "s#<FC>#${FC}#g" Makefile
sed -i "s#<FFLAGS>#${FFLAGS}#g" Makefile
sed -i "s#<HDF5_INCLUDE_DIR>#${HDF5_INCLUDE_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_DIR>#${HDF5_LIB_DIR}#g" Makefile
sed -i "s#<HDF5_LIB_FLAGS>#${HDF5_LIB_FLAGS}#g" Makefile
echo "Building MAPFL...."
make 1>build.log 2>build.err
echo "Copying MAPFL executable from SRC to BIN..."
cp ${MAPFL_HOME}/src/mapfl ${MAPFL_HOME}/bin/mapfl
echo "Done!"


export ROOTSYS='/opt/root/5.34-36-SL7-64-gcc6.1.0'

export TOP_DIR=`pwd | sed -e 's/\/nfs//'`
export extern=/data/sleipnir/PANDA/PWA
export LD_LIBRARY_PATH=$TOP_DIR/lib:$extern/libgcc610:$ROOTSYS/lib
export BOOST_BUILD_PATH=/data/sleipnir/PANDA/boost-build
export PATH=$BOOST_BUILD_PATH/bin:$TOP_DIR/bin:$ROOTSYS/bin:$PATH
export JULIA_DIR=/usr/local/julia
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$JULIA_DIR/lib/julia:$JULIA_DIR/lib
export JULIA_HOME=$JULIA_DIR/bin

rm Jamroot
ln -s JamrootSL7 Jamroot
rm Event/Jamfile
ln -s JamfilePython27 Event/Jamfile

echo "ROOTSYS: " $ROOTSYS 
echo "BOOST_BUILD_PATH: " $BOOST_BUILD_PATH  

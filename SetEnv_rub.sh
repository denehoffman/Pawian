source /opt/root/6-14.02-SL7-64-gcc9.2.0/bin/thisroot.sh

export TOP_DIR=`pwd | sed -e 's/\/nfs//'`
export extern=/data/sleipnir/PANDA/PWA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TOP_DIR/lib:$extern/libgcc920-opt
export BOOST_BUILD_PATH=/data/sleipnir/PANDA/boost-build
export PATH=$BOOST_BUILD_PATH/bin:$TOP_DIR/bin:$PATH

rm Jamroot
ln -s JamrootSL7 Jamroot
rm Event/Jamfile
ln -s JamfilePython27 Event/Jamfile

echo "ROOTSYS: " $ROOTSYS 
echo "BOOST_BUILD_PATH: " $BOOST_BUILD_PATH  

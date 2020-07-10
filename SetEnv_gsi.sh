export ROOTSYS=/lustre/nyx/panda/malte/extern/rootInstall

export TOP_DIR=`pwd | sed -e 's/\/nfs//'`
export extern=/lustre/nyx/panda/malte/extern
export LD_LIBRARY_PATH=$TOP_DIR/lib:$extern/lib:$ROOTSYS/lib
export BOOST_BUILD_PATH=/lustre/nyx/panda/malte/extern/boost-build/
export PATH=$ROOTSYS/bin:$TOP_DIR/bin:$PATH

rm Jamroot
ln -s JamrootSL6_gsi Jamroot

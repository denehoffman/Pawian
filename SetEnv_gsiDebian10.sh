#export ROOTSYS=/lustre/panda/mkuessne/extern/root-6.24.04/obj
#export ROOTSYS=/lustre/nyx/panda/malte/extern/rootInstall
export ROOTSYS=/lustre/panda/bkopf/extern/rootInstall/root/obj
export TOP_DIR=`pwd | sed -e 's/\/nfs//'`
export extern=/lustre/panda/bkopf/extern/
export LD_LIBRARY_PATH=$TOP_DIR/lib:$extern/lib:$ROOTSYS/lib
export BOOST_BUILD_PATH=/lustre/panda/bkopf/extern/boost_1_61_0
export PATH=$ROOTSYS/bin:$TOP_DIR/bin:$PATH

rm Jamroot
ln -s JamrootDebian10_gsi Jamroot

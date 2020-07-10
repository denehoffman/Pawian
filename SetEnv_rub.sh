source /opt/root/6-18.04-Centos8-gcc9.2.0/bin/thisroot.sh

export TOP_DIR=`pwd | sed -e 's/\/nfs//'`
export extern=/data/sleipnir/PANDA/PWA
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TOP_DIR/lib:$extern/libgcc920-opt
export BOOST_BUILD_PATH=/usr/share/boost-build
export PATH=$TOP_DIR/bin:$PATH

rm Jamroot
ln -s JamrootCentos8 Jamroot

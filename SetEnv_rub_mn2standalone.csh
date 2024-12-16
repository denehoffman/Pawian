# source /opt/root/6-28.04-AL9.2-gcc12.2.0/bin/thisroot.csh
source /opt/root/6-32.04-AL9.4-gcc12.2.0/bin/thisroot.csh
setenv BOOST_BUILD_PATH /usr/share/boost-build
setenv MINUIT2_STANDALONE /data/iltschi/pawian_utils/Minuit2/Minuit2-6.33.01_AlmaLinux9.4_gcc12.2.0
setenv KMAT_DIR /data/duldul/bertram/KMatStore/
setenv EVT_DIR /data/duldul/bertram/EvtStore/
setenv TOP_DIR `pwd | sed -e 's/\/nfs//'`
setenv extern /data/sleipnir/PANDA/PWA
setenv LD_LIBRARY_PATH ${MINUIT2_STANDALONE}/lib:${TOP_DIR}/lib:${LD_LIBRARY_PATH}
setenv PATH ${TOP_DIR}/bin:${PATH}

rm Jamroot
ln -s Jamroot_rub_AL9_mn2standalone Jamroot

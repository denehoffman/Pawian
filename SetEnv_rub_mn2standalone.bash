#### check for existing pawian environment
NEWTOPDIR=$(realpath "${BASH_SOURCE[0]}")
NEWTOPDIR="${NEWTOPDIR%/*}"
if [[ -z "$TOP_DIR" ]]; then
    TOP_DIR="${NEWTOPDIR}"
else
    echo "pawian env already set to ${TOP_DIR}"
    read -p "setup new env to ${NEWTOPDIR}? (yY/nN): " choice
    if [[ "$choice" == "y" ]] || [[ "$choice" == "Y" ]]; then
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH/${MINUIT2_STANDALONE}\/lib:${TOP_DIR}\/lib/}
        LD_LIBRARY_PATH=${LD_LIBRARY_PATH/::/:}
        if [[ "${LD_LIBRARY_PATH:0:1}" == ":" ]]; then
            LD_LIBRARY_PATH=${LD_LIBRARY_PATH:1}
        fi
        PATH=${PATH/${TOP_DIR}\/bin/}
        PATH=${PATH/::/:}
        if [[ "${PATH: -1}" == ":" ]]; then
            PATH=${PATH:0: -1}
        fi
        TOP_DIR="$NEWTOPDIR"
    elif [[ "$choice" == "n" ]] || [[ "$choice" == "N" ]]; then
        echo "NO, pawian env remains to $TOP_DIR"
        return
    else
        echo "invalid input, exiting"
        return
    fi
fi

############################
#### user defined variables
############################
# path to root installation (will only be used if root is not already sourced -> ROOTSYS not set)
ROOT_DIR="/opt/root/6-28.04-AL9.2-gcc12.2.0"
# path to Minuit2 package
MINUIT2_STANDALONE="/data/iltschi/pawian_utils/Minuit2/Minuit2-6.33.01_AlmaLinux9.4_gcc12.2.0"
# path to boost build
BOOST_BUILD_PATH="/usr/share/boost-build"
# pawian internal (path to KMatrix store)
KMAT_DIR="/data/duldul/bertram/KMatStore/"
# pawian internal (path to Event store)
EVT_DIR="/data/duldul/bertram/EvtStore/"
# set Jamfile for bjam/b2 (a symlink Jamroot -> $JAM_FILE will be created)
JAM_FILE="Jamroot_rub_AL9_mn2standalone"

#### source root
source ${ROOT_DIR}/bin/thisroot.sh

#### prepend minuit2 and pawian libs to LD_LIBRARY_PATH
LD_LIBRARY_PATH="${MINUIT2_STANDALONE}/lib:${TOP_DIR}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

#### append PATH
PATH="${PATH:+${PATH}:}${TOP_DIR}/bin"

#### export variables
export TOP_DIR
export KMAT_DIR
export EVT_DIR
export BOOST_BUILD_PATH
export MINUIT2_STANDALONE
export LD_LIBRARY_PATH
export PATH

# create symlink to Jamroot file (remove old link before)
if [[ -L "${TOP_DIR}/Jamroot" ]]; then
    jam_link_target=$(readlink "${TOP_DIR}/Jamroot")
    if [[ "${jam_link_target}" != "${JAM_FILE}" ]]; then
        rm ${TOP_DIR}/Jamroot
        ln -s ${JAM_FILE} ${TOP_DIR}/Jamroot
    fi
else
    ln -s "${JAM_FILE}" ${TOP_DIR}/Jamroot
fi

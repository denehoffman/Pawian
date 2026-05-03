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
ROOT_DIR="/home/dene/opt/root/root_install"
# path to Minuit2 package
MINUIT2_STANDALONE="$HOME/.local"
# path to boost build
# BOOST_BUILD_PATH="/usr/share/boost-build"
BOOST_BUILD_PATH="$HOME/opt/boost/1.70/tools/build"
# pawian internal (path to KMatrix store)
KMAT_DIR="/data/duldul/bertram/KMatStore/"
# pawian internal (path to Event store)
EVT_DIR="/data/duldul/bertram/EvtStore/"
# set Jamfile for bjam/b2 (a symlink Jamroot -> $JAM_FILE will be created)
JAM_FILE="Jamroot_andromeda"

#### source root
source ${ROOT_DIR}/bin/thisroot.sh

#### prepend minuit2 and pawian libs to LD_LIBRARY_PATH
LD_LIBRARY_PATH="${MINUIT2_STANDALONE}/lib:${TOP_DIR}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

LD_LIBRARY_PATH="$HOME/opt/boost/1.70/lib:$LD_LIBRARY_PATH"
CPLUS_INCLUDE_PATH="$HOME/opt/boost/1.70/include:$CPLUS_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH

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

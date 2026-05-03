#### check for existing pawian environment
NEWTOPDIR="$(realpath "${(%):-%x}")"
NEWTOPDIR="${NEWTOPDIR:h}"

if [[ -z "$TOP_DIR" ]]; then
    TOP_DIR="$NEWTOPDIR"
else
    echo "pawian env already set to ${TOP_DIR}"
    echo -n "setup new env to ${NEWTOPDIR}? (yY/nN): "
    read choice

    if [[ "$choice" == "y" || "$choice" == "Y" ]]; then
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH/${MINUIT2_STANDALONE}\/lib:${TOP_DIR}\/lib/}"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH//::/:}"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH#:}"
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH%:}"

        PATH="${PATH/${TOP_DIR}\/bin/}"
        PATH="${PATH//::/:}"
        PATH="${PATH#:}"
        PATH="${PATH%:}"

        TOP_DIR="$NEWTOPDIR"
    elif [[ "$choice" == "n" || "$choice" == "N" ]]; then
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
ROOT_DIR="/home/dene/opt/root/root_install"
MINUIT2_STANDALONE="$HOME/.local"
BOOST_BUILD_PATH="$HOME/opt/boost/1.70/tools/build"
KMAT_DIR="/data/duldul/bertram/KMatStore/"
EVT_DIR="/data/duldul/bertram/EvtStore/"
JAM_FILE="Jamroot_andromeda"

#### source root
source "${ROOT_DIR}/bin/thisroot.sh"

#### prepend boost, minuit2, and pawian libs
LD_LIBRARY_PATH="$HOME/opt/boost/1.70/lib:${MINUIT2_STANDALONE}/lib:${TOP_DIR}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CPLUS_INCLUDE_PATH="$HOME/opt/boost/1.70/include${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}"

#### append PATH
PATH="${PATH:+${PATH}:}${TOP_DIR}/bin"

#### export variables
export TOP_DIR
export KMAT_DIR
export EVT_DIR
export BOOST_BUILD_PATH
export MINUIT2_STANDALONE
export LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH
export PATH

#### create Jamroot symlink
if [[ -L "${TOP_DIR}/Jamroot" ]]; then
    jam_link_target="$(readlink "${TOP_DIR}/Jamroot")"
    if [[ "$jam_link_target" != "$JAM_FILE" ]]; then
        rm "${TOP_DIR}/Jamroot"
        ln -s "$JAM_FILE" "${TOP_DIR}/Jamroot"
    fi
else
    ln -s "$JAM_FILE" "${TOP_DIR}/Jamroot"
fi

#
# source this script to add concoord to the environment, like:
#
# source CONCOORDRC.bash
#
#this is the main directory, edit this line:
#
export CONCOORDDIR=$HOME/concoord2.1

#set the environment vaiables
export CONCOORDBIN=$CONCOORDDIR/bin
export CONCOORDLIB=$CONCOORDDIR/lib

# add binary path to $path
export PATH="$CONCOORDBIN":"$PATH"

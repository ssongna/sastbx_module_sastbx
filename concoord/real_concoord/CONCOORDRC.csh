#
# source this script to add concoord to the environment, like:
#
# source CONCOORDRC.csh  
#
#this is the main directory, edit this line:
#
setenv CONCOORDDIR $HOME/programs/concoord2.1

#set the environment vaiables
setenv CONCOORDBIN   $CONCOORDDIR/bin
setenv CONCOORDLIB   $CONCOORDDIR/lib

# add binary path to $path
setenv  PATH    "$CONCOORDBIN":"$PATH"

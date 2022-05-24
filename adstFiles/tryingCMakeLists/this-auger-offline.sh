export AUGEROFFLINEROOT=$(cd $(dirname $BASH_SOURCE)/.. ; pwd)
echo "export AUGEROFFLINEROOT=$AUGEROFFLINEROOT;"
if [ -e $AUGEROFFLINEROOT/bin/auger-offline-config ] ; then
  $AUGEROFFLINEROOT/bin/auger-offline-config --env-sh
  eval $($AUGEROFFLINEROOT/bin/auger-offline-config --env-sh)
  export PATH=$ROOTSYS/bin:$PATH
else
  echo "Cannot find auger-offline-config!"
  return 1
fi

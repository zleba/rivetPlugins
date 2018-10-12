#!/bin/bash
address=$PWD/../
#/nfs/dust/cms/user/zlebcr/powheg/armando/pythia/pythiaLocal/share/Pythia8/examples
cd $TMP
echo $PWD


name=bjetsNew

source $address/setup49.sh
fifo="hepmc${1}.fifo"
mkfifo $fifo

cp $address/Rivet${name}.so .
cp  $address/main41 .
cp  $address/main41.cmnd .

./main41 $fifo $1  &
rivet --pwd -a $name $fifo

cp Rivet.yoda  $address/farm/histos/Rivet${1}.yoda
rm $fifo

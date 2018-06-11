#!/bin/bash
address=$PWD/../
#/nfs/dust/cms/user/zlebcr/powheg/armando/pythia/pythiaLocal/share/Pythia8/examples
cd $TMP
echo $PWD


name=bjetsDeltaHL

source $address/setup62.sh
fifo="hepmc${1}.fifo"
mkfifo $fifo

cp $address/Rivet${name}.so .
cp  $address/main41 .

./main41 $fifo $1  &
rivet --pwd -a $name $fifo

cp Rivet.yoda  $address/farm/histos/Rivet${1}.yoda
rm $fifo

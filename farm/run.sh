#!/bin/bash
address=/nfs/dust/cms/user/zlebcr/powheg/armando/pythia/pythiaLocal/share/Pythia8/examples
cd $TMP
echo $PWD


source /nfs/dust/cms/user/zlebcr/powheg/armando/setup.sh
fifo="hepmc${1}.fifo"
mkfifo $fifo

$address/main41 $fifo $1  &
ln -s $address/RivetBjetsHL.so
rivet --pwd -a CMS_2011_S9086218_BJetak413TeV $fifo

cp Rivet.yoda  $address/farmnew/histos/Rivet${1}.yoda
rm $fifo

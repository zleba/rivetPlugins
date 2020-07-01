#!/bin/bash

mkdir -p /tmp/zlebcr/
fifo="/tmp/zlebcr/hepmc.fifo"
rm -f  $fifo
mkfifo $fifo

./main41 test.hepmc 10 #&
#./main41 $fifo 10 &
#rivet --pwd -a  ATLAS_2016_CONF_2016_092  $fifo

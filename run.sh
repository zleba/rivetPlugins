#!/bin/bash

mkdir -p /tmp/zlebcr/
fifo="/tmp/zlebcr/hepmc.fifo"
rm -f  $fifo
mkfifo $fifo

./main41 $fifo 10 &
rivet --pwd -a  bjetsNew $fifo

#!/bin/bash


fifo="hepmc.fifo"
mkfifo $fifo

./main41 $fifo 10 &
rivet --pwd -a CMS_2011_S9086218_BJetak413TeV $fifo

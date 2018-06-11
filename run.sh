#!/bin/bash


fifo="hepmc.fifo"
mkfifo $fifo

./main41 $fifo 10 &
rivet --pwd -a  bjetsDeltaHL $fifo

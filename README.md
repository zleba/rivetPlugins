# rivetPlugins

Simple setup for machine with cvmfs mounted (e.g. DESY naf or CERN lxplus):
```
git clone git@github.com:zleba/rivetPlugins.git
cd rivetPlugins
. ./setup62.sh
make main41 RivetbjetsDeltaHL.so
```

Then run the test localy using:
```
./run.sh
```
Or submit to farm using the submit file in the farm directory:
```
cd farm
condor_submit jobs.submit
```

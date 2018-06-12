


pythLib=/cvmfs/sft.cern.ch/lcg/releases/LCG_93a/MCGenerators/pythia8/230/x86_64-slc6-gcc62-opt/lib/
pythInc=/cvmfs/sft.cern.ch/lcg/releases/LCG_93a/MCGenerators/pythia8/230/x86_64-slc6-gcc62-opt/include/

hepMCinc=/cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc62-opt/include/
hepMClib=/cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc62-opt/lib


main41: main41.cc
	g++ $< -o $@ -I${hepMCinc} -I${pythInc} -O2  -pedantic -W -Wall -Wshadow -fPIC -L${pythLib} -Wl,-rpath,${pythLib} -lpythia8 -ldl \
     -L${hepMClib} -Wl,-rpath,${hepMClib} -lHepMC


RivetbjetsHL.so: bjetsHL.cc
	rivet-buildplugin $@ $<

RivetbjetsDeltaHL.so: bjetsDeltaHL.cc
	rivet-buildplugin $@ $<

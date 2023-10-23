#NuSTEC FSI hackanalysis

Requires NuHepMC: Build like:

```bash
git clone git@github.com:NuHepMC/cpputils.git
cd cpputils; mkdir build; cd build;

#builtin is important currently (23/10/23) as we are tracking 
#some new upstream features not in a release
cmake .. -DBUILTIN_HEPMC3=ON -DCMAKE_INSTALL_PREFIX=$(readlink -f Linux)
make install -j 8

#set up the NuHepMC environment
eval $(Linux/bin/NuHepMC-config --env)

#now we can build the tools in this directory
cd ../../

#this --build flag just invokes g++ with some useful flags for linking straight to NuHepMC tools
# anything after the analysis name is forwarded straight to the compiler
#we need root and if HepMC3 picked up the compression libs, then we need to pass those DSOs on the CLI
NuHepMC-config --build nustecana.cxx $(root-config --glibs --cflags) -llzma -lz -lbz2 -g -O0 -lfmt
NuHepMC-config --build dumptopy.cxx $(root-config --glibs --cflags) -llzma -lz -lbz2 -g -O0 -lfmt

#run the analysis
./nustecana <inp.hepmc3> <outputfile.root>
#turn the root files into an eval-able python literal that numpy can parse nicely
./dumptopy <outputfile.root> <generator tag> > hists.pynp
```

See an example [output](./hists.pynp).
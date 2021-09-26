# HF-Motecarlo Simulation
The codes of MonteCarlo method which aims to give the simulation of MuSEUM experiment; particle/nuclear physics experiment to measure the Hyperfine Splitting of Muonium at J-PARC. ROOT and GEANT4 is used here. <br>
GEANT4 is the simulation package of radiation. It is used in the fields of physics, etc. medical...
https://geant4.web.cern.ch/node/1<br>
ROOT is the analysis software flamework for high energy physis experiment. It is provided by cern.
https://root.cern/<br>
https://g-2.kek.jp/gakusai/research_MuHFS.html

1. RF.cc and RF.hh.<br>
To calculate the TMmode of microware in the cavity. We use TM110 and TM210 mode to demo the RF respectively.

2. magnet.cc and magnet.hh.<br>
To calculate the magnetic field about 1.7 T of superconductive magnet. This magnet will use at MuSEUM experiment and muon(g-2) experiment in future.

3. stop.cc and stop.hh.<br>
To calculate the muon stopping distribution in the RFcavity Using Geant4. 

4. simulator.cc and simulator.hh.<br>
To calculate the effective microwave distribution by superposition of RF field and muons distribution.

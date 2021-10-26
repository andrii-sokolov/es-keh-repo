# Modelling of Frequency Up-Conversion on a basis of Gap-Closing Electrostatic KEH 

This software is a part of a PhD project "Investigation and Optimisation of Kinetic Energy Harvesters with Nonlinear Electromechanical Coupling Mechanisms". It contains the procedure of the reconstruction on the resonance curves reconstruction and analysis of the waveforms. 

## Requirements

1. C++ compiling:
  * gcc >= 7.8
  * make >= 4.3
  * boost >= 1.65
2. Python:
  * python >= 3.7
  * jupyter lab >= 3.0
  * numpy >= 1.18

## Usage

Compilation of the C++ software: 
```
make
```
Execution of the C++ software: 
```
make execute
```
Plot the graphs: 
```
make plot
```
Clear all temporary files:
```
make clear
```
All analysis software is in the `pysrc` folder

#!/bin/bash -i
rm qrot-plot*.png 
clang++ -o qrot -O2 -Wall qrot.cpp 
./qrot --dt=.1e-5 --steps=25e6 --omega=10,0.01,0 --verbosity=1 --trace-file=axis1.csv --trace-step=1000
./qrot --dt=.1e-5 --steps=25e6 --omega=0,10,0.01 --verbosity=1 --trace-file=axis2.csv --trace-step=1000
./qrot --dt=.1e-5 --steps=25e6 --omega=0.01,0,10 --verbosity=1 --trace-file=axis3.csv --trace-step=1000
python qrot-plot.py --trace-files axis1.csv axis2.csv axis3.csv --figext png

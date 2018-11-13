#! /bin/bash
set -o nounset

DATADIR=$HOME/data/HSC/brightStarMasks

HOST=$( hostname -s )

TMPDIR=$HOME/data/tmp
STILTS='java -Xmx8192M -jar /Users/coupon/local/bin/stilts.jar'

# default mpirun and swot locations
MPIRUN=mpirun
SWOT=$HOME/local/source/GitHub/swot/bin/swot
VENICE=$HOME/local/source/GitHub/venice/bin/venice

# default number of cores
NP=2
if [[ "$HOST" =~ "ojingo" ]]; then
   # compile with make MPICC=/opt/openmpi-1.8.6_clang/bin/mpicc
   MPIRUN=/opt/openmpi-1.8.6_clang/bin/mpirun
   NP=8
fi
if [[ "$HOST" =~ "hayabusa" ]]; then
   MPIRUN=/opt/openmpi-1.8.5/bin/mpirun
   NP=4
fi

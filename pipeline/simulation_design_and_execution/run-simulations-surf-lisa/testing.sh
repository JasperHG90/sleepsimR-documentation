#!/bin/bash

TMPDIR=~/tmp
mkdir $TMPDIR

# Copy singularity container / R script to node
cp -r $HOME/input_dir "$TMPDIR"/input_dir

# Make the output directory for storing entire models
mkdir "$TMPDIR"/var && mkdir "$TMPDIR"/var/sleepsimR

#Determining the number of processors in the system
NPROC=6

# CD into the tempdir
cd $TMPDIR/input_dir
# Run at most 1.000 * 16 == 16.000 iterations
for i in `seq 1 2`; do
  python3 singularity-run-simulation.py 'http://35.214.134.76/6532a6c0213d9ac79f2101838e40e041' 'e5bafe0847c32243e909197c20f6ab44' '46104911225654f814159cd940e4364b' 'singularity run -B ~/scratch/var/sleepsimR:/var/sleepsimR sleepsimR-run-v03-beta.simg' --iterations=$NPROC && rsync -avhu --progress "$TMPDIR"/var/sleepsimR/ $HOME/sleepsimrmodels/
done
wait

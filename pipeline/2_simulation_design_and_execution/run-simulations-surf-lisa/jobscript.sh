#!/bin/bash

# Adapted from:
#   - example: Single node, parallel execution of a serial program on all cores - identical input, using scratch
#   - url: https://userinfo.surfsara.nl/systems/lisa/user-guide/creating-and-running-jobs#example-job-scripts
#   - from: surfsara

# Use 16 cores on a single node
#SBATCH -N 1

# Wall clock time. Node is allocated for 1 hr
#SBATCH -t 30:00:00

# Request normal partition
#SBATCH -p normal

# Mail at beginning and end of job
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jasperginn@gmail.com

# Copy singularity container / R script to node
cp -r $HOME/input_dir "$TMPDIR"/input_dir

# Make the output directory for storing entire models
mkdir "$TMPDIR"/var && mkdir "$TMPDIR"/var/sleepsimR

#Determining the number of processors in the system
NPROC=`nproc --all`

# Run job. We are on a single node with 16 cores, so we want to push out 16 jobs at once.
# ... Make another for-loop here to run multiple programs
# We are taking three actions:
#  (1) CD into the temporary directory with singularity image and application
#      unlike docker, singularity does not store the program in its layers. 
#      It binds the current directory to the container.
#  (2) Run the singularity program. Attach another folder (/var/sleepsimR) that is
#      used to store full models.
#  (3) When program is done. Sync the /var/sleepsimR folder with the folder $HOME/sleepsimrmodels
#      to ensure that any models that are stored (only happens in 5% of all iterations) are then
#      also stored persistently on disk. Two other reasons for rsync:
#       (a) We don;t know if a model file exists, and we don't want to use a conditional statement here
#       (b) We don't want to wait until the end of the program. We have at most 120hrs of computing time.
#           However, earlier iterations of the program take much less time than later iterations. (it is
#           the difference between 10 subj. and e.g. 400 timesteps and 80 subjects and 3.200 timesteps).
#           So we are running this program in a double loop that basically ends when the wall time ends.
#           This means that the current iterations that are running when time is up will never be completed.
#           However, this is not really an issue because they are sent back to the pool if they don't finish
#           within a certain timespan.

# CD into the tempdir
cd $TMPDIR/input_dir
# Run at most 1.000 * 16 == 16.000 iterations
for i in `seq 1 1000`; do
  python3 singularity-run-simulation.py 'http://35.214.134.76/6532a6c0213d9ac79f2101838e40e041' 'e5bafe0847c32243e909197c20f6ab44' '46104911225654f814159cd940e4364b' 'singularity run -B /scratch/var/sleepsimR:/var/sleepsimR sleepsimR-run-v03-beta.simg' --iterations=$NPROC && rsync -avhu --progress "$TMPDIR"/var/sleepsimR/ $HOME/sleepsimrmodels/
done
wait

# Copy output back
# (Make conditional .. only if model has been saved!)
# Use rsync to sync these directories. 
# See https://unix.stackexchange.com/questions/149965/how-to-copy-merge-two-directories
# rsync -avhu --progress "$TMPDIR"/var/sleepsimR/ $HOME/sleepsimrmodels/

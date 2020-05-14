# Creating Slurm batch jobs on the LISA system

The LISA cluster uses Slurm to process batch jobs. Documentation can be found [here](https://userinfo.surfsara.nl/systems/lisa/user-guide/creating-and-running-jobs#example-job-scripts).

The job script that I use is enclosed in this folder ("jobscript.sh"). 

Open this file in your favorite text editor, and change the IP address to that of the "sleepsimR-api" program. Also change the password and username if you changed it. Finally, ensure that you are using the right name of the sleepsimR-run Singularity image. 

```shell
...
for i in `seq 1 1000`; do
  python3 singularity-run-simulation.py 'http://<YOUR-IP-ADDRESS>/6532a6c0213d9ac79f2101838e40e041' <USERNAME> <PASSWORD> 'singularity run -B /scratch/var/sleepsimR:/var/sleepsimR <SLEEPSIMR-RUN-IMAGE>.simg' --iterations=$NPROC && rsync -avhu --progress "$TMPDIR"/var/sleepsimR/ $HOME/sleepsimrmodels/
done
...
```

For example:

```shell
...
for i in `seq 1 1000`; do
  python3 singularity-run-simulation.py 'http://35.214.134.76/6532a6c0213d9ac79f2101838e40e041' 'e5bafe0847c32243e909197c20f6ab44' '46104911225654f814159cd940e4364b' 'singularity run -B /scratch/var/sleepsimR:/var/sleepsimR sleepsimR-run-version-1.3.simg' --iterations=$NPROC && rsync -avhu --progress "$TMPDIR"/var/sleepsimR/ $HOME/sleepsimrmodels/
done
...
```

Save the file and put it into the "input_dir" that also contains the Singularity image. Make sure to also copy the "singularity-run-simulation.py" file into the "inpur_dir" directory.

Next, SCP the entire "input_dir" folder to your LISA home folder.

Finally, boot into your LISA home directory and create the following folder:

```shell
mkdir sleepsimrmodels
```

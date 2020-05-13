# Running simulations on the LISA cluster

If you're on a fresh Ubuntu 18.04 OS, install the following prerequisites

```shell
sudo apt-get install -y python libarchive-dev
```

Download and install singularity 2.6 by following [the instructions](https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps)

Build the singularity file

```shell
singularity build sleepsimR-run-v03-beta.simg docker://jhginn/sleepsimr-run:version-1.2
```

Make folder sleepsimrmodels

```shell
mkdir ~/sleepsimrmodels
```

Make folder input_dir

```shell
mkdir ~/input_dir
```

clone repoistory 'sleepsimR-documentation'

```shell
git clone git@github.com:JasperHG90/sleepsimR-documentation.git
```

Clone the 'sleepsimR-run' repository

```shell
git clone git@github.com:JasperHG90/sleepsimR-run.git
```

Copy the app to the input dir folder

```shell
cp -r sleepsimR-run/app input_dir
```

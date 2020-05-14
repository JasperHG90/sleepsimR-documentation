# Create a Singularity image

First, download and install singularity 2.6 **on your own device** by following [the instructions](https://sylabs.io/guides/2.6/user-guide/quick_start.html#quick-installation-steps).

Next, download and unzip the version of "sleepsimR-run" that pairs with the API version you are running (see entries in table 1 in main README file):

```shell
wget https://zenodo.org/record/3727710/files/JasperHG90/sleepsimR-run-v1.3.zip?download=1 -O sleepsimR-run-v1.3.zip && unzip sleepsimR-run-v1.3.zip
```

Go into the directory and build the dockerfile:

```shell
cd JasperHG90-sleepsimR-run-198260a && docker build . -t jhginn/sleepsimr-run:version-1.3
```

On your device, set up a local docker registry:

```shell
# Start a docker registry
docker run -d -p 5000:5000 --restart=always --name registry registry:2
# Push local docker container to it
docker tag jhginn/sleepsimr-run:version-1.3 localhost:5000/jhginn/sleepsimr-run:version-1.3
docker push localhost:5000/jhginn/sleepsimr-run:version-1.3
```

Next, create a Singularity file by executing 

```shell
nano Singularity 
```

And copy-paste the following text into it:

```text
Bootstrap: docker
Registry: http://localhost:5000
Namespace:
From: jhginn/sleepsimr-run:version-1.3
```

Finally, build the Singularity image:

```shell
sudo SINGULARITY_NOHTTPS=1 singularity build sleepsimR-run-v1.3.simg Singularity
```

Stop the local registry

```shell
docker stop registry && docker rm registry
```

Move all required files into a folder called "input_dir":

```shell
mkdir input_dir && mv sleepsimR-run-v1.3.simg input_dir && mv JasperHG90-sleepsimR-run-198260a/app input_dir
```

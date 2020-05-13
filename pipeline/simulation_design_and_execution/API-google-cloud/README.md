# Setting up the sleepsimR-api program

Boot into the virtual machine and follow the instructions given [here](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04) to install docker and the instructions given [here](https://www.digitalocean.com/community/tutorials/how-to-install-docker-compose-on-ubuntu-18-04) to install docker compose. 

Next, download a version of the API from one of the Zenodo repositories and unzip it:

```shell
wget https://zenodo.org/record/3727773/files/JasperHG90/sleepsimR-api-v1.3.1.zip?download=1 -O sleepsimR-api-v1.3.1.zip && unzip sleepsimR-api-v1.3.1.zip
```

Build the docker application as follows (NB: change this if you use a different version):

```shell
cd JasperHG90-sleepsimR-api-878f033 && docker build . -t jhginn/sleepsimr-api:version-1.3.1 && cd
```

Create the following docker volumes:

```shell
docker volume create sleepsimr && docker volume create sleepsimrmodels
```

Next, download this research archive on the virtual machine:

```shell
git clone https://github.com/JasperHG90/sleepsimR-documentation.git
```

Copy the following folder to your home directory

```shell
cp -R ~/thesis-docs/pipeline/simulation_design_and_execution/API-google-cloud ~/docker-api
```

Enter this folder and open the "docker-compose.yml" file using the following command:

```shell
nano docker-compose.yml
```

Change that line that is enclosed in the red box in the image below to the version that you are using. For example, I change the line to `jhginn/sleepsimr-api:version-1.3.1`. 

Next, look for the lines that say "SLEEPSIMR_USR" and "SLEEPSIMR_PWD" and change these values. I strongly recommend that you do this. Hit \<ctrl\> + X and then \<ENTER\> to save your changes.

<figure>
  <img src="img/api-dc.png" alt="Pipeline" style="width:80%">
</figure> 

Now, execute the following command to start the API:

```shell
docker-compose up -d 
```

The API is now up and running on http://\<YOUR-VM-IP-ADDRESS\>/6532a6c0213d9ac79f2101838e40e041


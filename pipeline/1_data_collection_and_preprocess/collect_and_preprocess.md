# Data collection and preprocessing

First, download the sleepsimR-collect and sleepsimR-preprocess programs and unzip the contents in a directory of your choosing;

```shell
wget https://zenodo.org/record/3826081/files/JasperHG90/sleepsimR-collect-v0.3.zip?download=1 -O sleepsimR-collect-v0.3.zip && wget https://zenodo.org/record/3826082/files/JasperHG90/sleepsimR-preprocess-v0.3.zip?download=1 -O sleepsimR-preprocess-v0.3.zip && unzip sleepsimR-collect-v0.3.zip && unzip sleepsimR-preprocess-v0.3.zip
```

Enter each repository and build the docker programs:

```shell
cd JasperHG90-sleepsimR-collect-ea01015 && docker build . -t jhginn/sleepsimr-collect && cd .. && cd JasperHG90-sleepsimR-preprocess-02976b1 && docker build . -t jhginn/sleepsimr-preprocess & cd ..
```

Create a docker volume to store the data:

```shell
docker volume create sleepsimr_data
```

Then, download the first night of EDFX sleep data by executing the following command:

```shell
docker run --rm --mount source=sleepsimr_data,target=/var/sleepsimr_data jhginn/sleepsimr-collect 1:83 1 /var/sleepsimr_data/raw /var/sleepsimr_data/processed
```

Download the second night of EDFX sleep data by executing the following command:

```shell
docker run --rm --mount source=sleepsimr_data,target=/var/sleepsimr_data jhginn/sleepsimr-collect 1:83 2 /var/sleepsimr_data/raw /var/sleepsimr_data/processed
```

Depending on the speed of your internet connection, downloading this database can be quite slow. However, you can stop the program at any time and resume downloading the data using the above commands.

Next, we run the preprocessing container:

```shell
docker run --rm --mount source=sleepsimr_data,target=/var/sleepsimr_data jhginn/sleepsimr-preprocess /var/sleepsimr_data/processed /var/sleepsimr_data/final
```

To copy the data from the docker volume, execute the following commands:

```shell
docker run --mount source=sleepsimr_data,target=/var/sleepsimr_data --name helper busybox
```

```shell
docker cp helper:/var/sleepsimr_data/final/EEG_data_final.rds .
```

Then, remove the busybox container

```shell
docker stop helper && docker rm helper
```

The final preprocessing steps for the `EEG_data_final.rds` data file are conducted in the [sleepsimRdata R library](https://github.com/JasperHG90/sleepsimRdata/blob/master/data-raw/1_preprocess_sleep_data/preprocess_sleep_data.R).

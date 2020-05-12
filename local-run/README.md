# Running the simulations locally

Make sure the docker images are up-to-date

```shell
docker pull jhginn/sleepsimr-run:latest && docker pull jhginn/sleepsimr-api:latest
```

First, set up the API. This can be done by doing the following:

```shell
docker run --rm --mount source=sleepsimr,target=/var/sleepsimR -p 5002:5002 -e SLEEPSIMR_USER="myuser" -e SLEEPSIMR_PWD="mypassword" --name sleepsimr-api jhginn/sleepsimr-api:latest
```

Run a single iteration container by executing:

```shell
docker run --rm --mount source=sleepsimrmodels,target=/var/sleepsimR jhginn/sleepsimr-run:devel --username "myuser" --password "mypassword" --host "http://172.17.0.1:5002"
```

# sleepsimR-documentation

This repository contains the research archive accompanying my Master's thesis. 

## Thesis details

| Title | Sample Size Considerations for Bayesian Multilevel Hidden Markov Models: A Simulation Study and Application to Electroencephalogram (EEG) and Electrooculography (EOG) Data to Detect Sleep States |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| Author | Jasper Ginn |
| Supervisor | Dr. Emmeke Aarts |
| Department | Department of Methodology & Statistics, Utrecht University |
| Data collection (period) | November 2019 |
| Manuscript submitted (period) | May 2020 |

## Disclaimer

In the instructions given below, I assume that you have access to a bash-like terminal.

## Data management

### 1. Ethical approval

### 2. Summary of research paper

### 3. Project pipeline

### 4. Software requirements and installation

Docker / docker-compose
Singularity


### 5. Replicating the analysis

#### 5.1 Data collection & preprocessing

#### 5.2 Simulation design & execution

#### 5.3 Analysis & empirical application

### 6. (Final) datasets

### 7. Analyses and plots

## Privacy

## Permission and access

# Docker container --> running iterations

1. docker volume create sleepsimrmodels
2. cd run-simulations-docker-compose && docker-compose up --scale iteration=3 -d

# Setting up the API server

docker run --mount source=sleepsimR,target=/var/sleepsimR --name helper busybox

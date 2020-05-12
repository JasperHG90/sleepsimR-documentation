# sleepsimR-documentation

This repository contains the research archive accompanying my Master's thesis. 

## Thesis details

<table>
<thead>
  <tr>
    <th>Field</th>
    <th>Value<br></th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>Title</td>
    <td>Sample Size Considerations for Bayesian Multilevel Hidden Markov<br>Models: A Simulation Study and Application to Electroencephalogram<br>(EEG) and Electrooculography (EOG) Data to Detect Sleep States</td>
  </tr>
  <tr>
    <td>Author</td>
    <td>Jasper Ginn</td>
  </tr>
  <tr>
    <td>Supervisor</td>
    <td>Dr. Emmeke Aarts</td>
  </tr>
  <tr>
    <td>Department</td>
    <td>Department of Methodology &amp; Statistics, Utrecht University</td>
  </tr>
  <tr>
    <td>Data collection (period)<br></td>
    <td>November 2019</td>
  </tr>
  <tr>
    <td>Manuscript submitted (period)</td>
    <td>May 2020</td>
  </tr>
</tbody>
</table>

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

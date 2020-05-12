# sleepsimR-documentation

This repository contains the research archive accompanying my Master's thesis. 

## Thesis details

<table class="tg">
<thead>
  <tr>
    <th class="tg-0pky">Title</th>
    <th class="tg-dvpl">Sample Size Considerations for Bayesian Multilevel Hidden Markov<br>Models: A Simulation Study and Application to Electroencephalogram<br>(EEG) and Electrooculography (EOG) Data to Detect Sleep States</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">Author</td>
    <td class="tg-dvpl">Jasper Ginn</td>
  </tr>
  <tr>
    <td class="tg-0pky">Supervisor</td>
    <td class="tg-dvpl">Dr. Emmeke Aarts</td>
  </tr>
  <tr>
    <td class="tg-0pky">Department</td>
    <td class="tg-dvpl">Department of Methodology &amp; Statistics, Utrecht University</td>
  </tr>
  <tr>
    <td class="tg-0pky">Data collection (period)<br></td>
    <td class="tg-dvpl">November 2019</td>
  </tr>
  <tr>
    <td class="tg-0pky">Manuscript submitted (period)</td>
    <td class="tg-dvpl">May 2020</td>
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

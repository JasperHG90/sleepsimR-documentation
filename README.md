# sleepsimR-documentation

This repository contains the research archive accompanying my Master's thesis. 

## Thesis details and summary

Spurred in part by the ever-growing number of sensors and web-based methods of collecting data, the use of Intensive Longitudinal Data (ILD) is becoming more common in the social and behavioral sciences [1]. The ILD collected in these fields are often hypothesized to be the result of latent states (e.g. behavior, emotions), and the promise of ILD lies in its ability to capture the dynamics of these states as they unfold in time. In particular, by collecting data for multiple subjects, researchers can observe how such dynamics differ between subjects. The Bayesian Multilevel Hidden Markov Model (mHMM) is a relatively novel model that is suited to model the ILD of this kind while taking into account heterogeneity between subjects. In my thesis, I conduct a Monte Carlo simulation study in which I vary (1) the number of subjects, (2) the number of occasions and (3) the between-subject variance and study the effect of varying those quantities on the parameter estimates obtained from the mHMM.

Based on previous studies in the multilevel modeling literature generally and for ILD models more specifically (e.g. [2], [3], [4]), I hypothesize that the subject sample size is the most important determinant of parameter estimate quality. This expectation is largely corroborated in the thesis. However, the occasion sample size is found to be important to adequately model the latent state transitions. I discuss how data characteristics influence parameter estimates and provide recommendations to researchers seeking to apply the mHMM to their own data. 

## Overview of the software developed in this thesis

The figure below

 <figure>
  <img src="img/pipeline.png" alt="Pipeline" style="width:80%">
  <figcaption>Figure 1: thesis project pipeline from data collection to analysis of the results. I subdivide the pipeline into three major components: (1) data collection and preprocessing, (2)  The software that was developed throughout the project </figcaption>
</figure> 

<table>
<thead>
  <tr>
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

In the instructions given below, I assume that you have access to a bash-like terminal. I also note that you cannot use Docker on [Windows 10 Home](https://docs.docker.com/docker-for-windows/install/), although there is [a workaround](https://medium.com/@mbyfieldcameron/docker-on-windows-10-home-edition-c186c538dff3) that you can use.

## Data management

### 1. Ethical approval

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

## References

# Docker container --> running iterations

1. docker volume create sleepsimrmodels
2. cd run-simulations-docker-compose && docker-compose up --scale iteration=3 -d

# Setting up the API server

docker run --mount source=sleepsimR,target=/var/sleepsimR --name helper busybox

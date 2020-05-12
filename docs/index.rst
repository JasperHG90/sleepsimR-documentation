<img src="img/pipeline.png" alt="data pipeline thesis" width="800"> 

<table class="tg">
<thead>
  <tr>
    <th class="tg-fymr">Name</th>
    <th class="tg-fymr">URL</th>
    <th class="tg-fymr">DOI<br></th>
    <th class="tg-fymr">Version</th>
    <th class="tg-fymr">Description</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">sleepsimR</td>
    <td class="tg-0pky">https://github.com/JasperHG90/sleepsimR</td>
    <td class="tg-0pky">10.5281/zenodo.3804726</td>
    <td class="tg-0pky">0.5</td>
    <td class="tg-0pky">R library that contains convenience functions to generate simulation scenarios, to run an mHMM, to simulate datasets and to obtain MAP estimates.</td>
  </tr>
  <tr>
    <td class="tg-0pky">sleepsimRdata</td>
    <td class="tg-0pky">https://github.com/JasperHG90/sleepsimRdata</td>
    <td class="tg-0pky">XXX-XXX-XXX</td>
    <td class="tg-0pky">0.3</td>
    <td class="tg-0pky">R library that contains scripts used to generate simulation scenarios and to preprocess the simulation results. This library contains all raw data, summary statistics, simulation scenarios and simulation results obtained in the course of the study.</td>
  </tr>
  <tr>
    <td class="tg-0pky">sleepsimReval</td>
    <td class="tg-0pky">https://github.com/JasperHG90/sleepsimReval</td>
    <td class="tg-0pky">10.5281/zenodo.3804835</td>
    <td class="tg-0pky">0.8</td>
    <td class="tg-0pky">R library that contains simulation evaluation metrics and utility functions to parse simulation results and assess convergence of MCMC chains.</td>
  </tr>
  <tr>
    <td class="tg-0pky" rowspan="7">sleepsimR-api</td>
    <td class="tg-0pky" rowspan="7">https://github.com/JasperHG90/sleepsimR-api</td>
    <td class="tg-0pky">10.5281/zenodo.3727709</td>
    <td class="tg-0pky">1.3.1</td>
    <td class="tg-0pky">Docker application. Resource manager that I use to organize my simulation study. This version contains the first 7.000 iterations (approx. 48 iterations / scenario).</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3731364</td>
    <td class="tg-0pky">1.3.2</td>
    <td class="tg-0pky">This version contains the remaining +- 29.000 iterations of the simulation study.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3747158</td>
    <td class="tg-0pky">1.3.3</td>
    <td class="tg-0pky">This version contains +- 834 iterations (approx. 6 iterations / scenario) that are used to assess convergence of a small percentage of the models that are estimated in the simulation study.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3784910</td>
    <td class="tg-0pky">1.5.1</td>
    <td class="tg-0pky">This version contains the simulation instructions for baseline scenarios 1 and 2.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3784934</td>
    <td class="tg-0pky">1.5.2</td>
    <td class="tg-0pky">This version contains the simulation instructions for baseline scenario 5.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3786154</td>
    <td class="tg-0pky">1.5.3</td>
    <td class="tg-0pky">This version contains the simulation instructions for baseline scenario 3.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3793259</td>
    <td class="tg-0pky">1.5.4</td>
    <td class="tg-0pky">This version contains the simulation instructions for baseline scenario 4.</td>
  </tr>
  <tr>
    <td class="tg-0pky">sleepsimRapiClient</td>
    <td class="tg-0pky">https://github.com/JasperHG90/sleepsimRapiClient</td>
    <td class="tg-0pky">10.5281/zenodo.3805052</td>
    <td class="tg-0pky">1.0</td>
    <td class="tg-0pky">R library that contains convenience functions to retrieve simulation parameters from the resource manager and allows you to send back the results.<br></td>
  </tr>
  <tr>
    <td class="tg-0pky" rowspan="3">sleepsimR-run</td>
    <td class="tg-0pky" rowspan="3">https://github.com/JasperHG90/sleepsimR-run</td>
    <td class="tg-0pky">10.5281/zenodo.3727710</td>
    <td class="tg-0pky">1.3</td>
    <td class="tg-0pky">Docker application. Contains the R script that I use to run a single iteration of the simulation study. This version is used with versions 1.3.1, 1.3.2, 1.3.3 of the sleepsimR-api program.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3778191</td>
    <td class="tg-0pky">1.5.1</td>
    <td class="tg-0pky">This version is used with version 1.5.1 (Scenario 1)  of the sleepsimR-api program. In this version, the emission distribution means are spread out.</td>
  </tr>
  <tr>
    <td class="tg-0pky">10.5281/zenodo.3778195</td>
    <td class="tg-0pky">1.5.2</td>
    <td class="tg-0pky">This version is used with versions 1.5.1 (Scenario 2), 1.5.3 and 1.5.4 of the sleepsimR-api program. In this version, the emission distribution means are spread out and the self-transition probabilities are lowered.</td>
  </tr>
  <tr>
    <td class="tg-0pky">sleepsimR-sleepdata-analysis</td>
    <td class="tg-0pky">https://github.com/JasperHG90/sleepsimR-sleepdata-analysis</td>
    <td class="tg-0pky">XXX.XXX.XXX</td>
    <td class="tg-0pky">1.0</td>
    <td class="tg-0pky">Docker application. This program allows the user to run a single chain of the model used in the empirical analysis.</td>
  </tr>
</tbody>
</table>

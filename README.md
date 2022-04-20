# Masters internship project code

## Introduction
Our project investigates if brain volume has a confounding effect on brain connectivity, and if so, how is it seen in the constructed graph metrics? 
The aim is to determine whether an individual’s brain volume can affect their brain connectivity and act as a covariate in functional connectivity studies.

## Research design 
First we compute graph metrics (clustering coefficients and path lengths) of subjects who are taken from the HCP dataset. 
Then we order the subjects according to their brain volumes. After this we investigate whether there is a trend seen in the computed graph metrics. 

## What's done
1. Data is cleaned, outliers removed
2. Clustering coefficient computation
3. Shortest path length computation
4. Random reference networks computation
5. Normalized clustering coefficient computation
6. Normalized shortest path length computation
7. Small-world index computation

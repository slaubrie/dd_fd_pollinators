# dd_fd_pollinators
Code for analyzing annual plant data. Data were collected in an experiment to determine how plant communities are regulated in Western Australia 2020 

# Title of Dataset
---

Patterns of frequency and density dependence are highly variable in diverse annual flowering plant communities

## abstract
Applications of ecological theory to natural communities often assume that competitive, negative density-dependent processes are the only type of interaction important for diversity maintenance. Recent advances suggest that positive interactions within trophic levels (e.g. plant-plant) may also affect plant coexistence. Though positive plant-plant interactions theoretically might result in positive or nonmonotonic frequency or density dependence (FD/DD), less is known about how commonly these patterns occur, or which ecological processes might result in such patterns in natural plant communities. In this study we test for signals of variable frequency and density dependence in annual flowering plant communities in Western Australia, and search for evidence that interactions among plants during flowering might induce positive or nonmonotonic FD/DD in flowering plants. Using four common annual wildflower species, we ask if plant fecundity exhibits positive or nonmonotonic FD/DD, and if pollinator-mediated plant-plant interactions during flowering change patterns of FD/DD relative to pollinator-independent plant interactions. Three species exhibited nonmonotonic (hump-shaped) density dependence, and only one species experienced strictly negative density dependence. Each species exhibited a different pattern of frequency dependence (positive, negative, weakly nonmonotonic, and no detectable frequency dependence). Pollinator-mediated plant-plant interactions during flowering induced both nonmonotonic density dependence and negative frequency dependence in one species. Importantly, the extent of variation in FD/DD observed in our study brings into question the dominance of negative density and frequency dependence in theory, suggesting instead that demographic responses of plants to their communities fall along a continuum of possible density- and frequency-dependent patterns.

## Description of the data and file structure

There are three csv data files. 
The main file is pfdata_2020.csv and is used to analyze density and frequency dependence. The corresponding code is called 'densities_analysis_GLMM_globalModelOnly.R' 
The file Bagged_Plants2020.csv is used to analyze pollen limitation in the four species. The corresponding code is stored in 'pollen_lim_2020.R'
The file 'ring_richnessData_2020.csv' contains neighborhood richness data around the focal plants in the study. The corresponding code is stored in 'nmds_exra.R'

## Sharing/Access information


Links to other publicly accessible locations of the data:
  * https://github.com/slaubrie/dd_fd_pollinators


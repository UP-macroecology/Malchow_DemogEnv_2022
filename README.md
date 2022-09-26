# Demography-environment relationships improve mechanistic understanding of range dynamics under climate change

This is supplementory material to the publication:

*"Demography-environment relationships improve mechanistic understanding of range dynamics under climate change"*

by Anne-Kathleen Malchow, Florian Hartig, Jette Reeg, Marc Kéry, and Damaris Zurell.

It contains the R scripts and data needed to reproduce all results. The code uses a developmental version of the R package *RangeShiftR*, which can be found in the RangeShiftR repo under the development tag [v.1.1-beta.0](https://github.com/RangeShifter/RangeShiftR-package/releases/tag/v.1.1-beta.0).


**Funding**

AM and DZ were supported by Deutsche Forschungsgemeinschaft (DFG) under grant agreement No. ZU 361/1-1.


**Data provision**

The study uses data from the Swiss breeding bird survey ([Schmid et al., 2004](#1)) and the Swiss breeding bird index ([Knaus et al., 2022](#2)) provided by the Swiss ornithological institute, Sempach.
The repository contains processed versions of this data.

<a id="1"></a>
Schmid, H., Zbinden, N., and Keller, V. (2004). Überwachung der Bestandsentwicklung häufiger Brutvögel in der Schweiz. Swiss Ornithological Institute.

<a id="2"></a>
Knaus, P., Strebel, N., & Sattler, T. (2022). The State of Birds in Switzerland 2022. Swiss Ornithological Institute. http://www.vogelwarte.ch/state

---

## Folder *data*

All input data required to run the model calibration: Raster maps of climate variables, the habitat models and derived maps, the initial distribution models and spatially aggregated counts from the Swiss breeding bird survey.


## Folder *scripts*

#### *scripts/data_prep*

Scripts to generate the input data from the original data sources, i.e. Chelsa climate v2.1 and CORINE land cover.

#### *scripts/calibration*

The scripts to run the model calibration. The routine is contained in *calibration_main.R*, which loads functions defined in the other scripts.

#### *scripts/analysis*

The scripts to process and analyse the MCMC chains produced in the calibration. Extracts the demography-climate relationships and performs prior/posterior predicitons. Creates all plots. All routines are contained in "analysis_main.R*, which loads functions defined in the other scripts.

## Folder *model*

Contains the RangeShiftR folder structure.

---

## Session info

The calculations were run under the following setup:

**R version 4.0.4 (2021-02-15)**  
Platform: x86_64-pc-linux-gnu (64-bit)  
Running under: Debian GNU/Linux 11 (bullseye)  

**Matrix products:**  
default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

**Locale:**  
 LC_CTYPE=en_US.UTF-8; LC_NUMERIC=C; LC_TIME=C; LC_COLLATE=en_US.UTF-8; LC_MONETARY=C; LC_MESSAGES=en_US.UTF-8; LC_PAPER=C; LC_NAME=C; LC_ADDRESS=C; LC_TELEPHONE=C; LC_MEASUREMENT=C; LC_IDENTIFICATION=C  

**Attached base packages:**  
 stats; graphics; grDevices; utils; datasets; methods; base  

**Other attached packages:**  
 RangeShiftR_1.0.3; raster_3.5-21; sp_1.5-0  

**Loaded via a namespace (and not attached):**  
 compiler_4.0.4; rgdal_1.5-27; tools_4.0.4; Rcpp_1.0.7; codetools_0.2-18; grid_4.0.4; rbibutils_2.2.8; Rdpack_2.3.1; lattice_0.20-41; terra_1.5-34 


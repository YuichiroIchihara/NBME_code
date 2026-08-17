# Code accompanying the manuscript

This repository contains the R scripts used for the analyses reported in the manuscript.

## Contents

- BBBscore: statistical analysis of BBB scores
- bulkRNAseq: bulk RNA-seq analysis
- demo: example dataset and demonstration analysis
- Gastrocnemius_weight: statistical analysis of gastrocnemius muscle weight
- Histology: statistical analysis of histological data
- Kinematics: kinematic and PCA analyses
- MEP: statistical analysis of motor-evoked potential data
- snRNAseq: single-nucleus RNA-seq analysis

## Software requirements

Analyses were performed using R version 4.4.2 in RStudio on macOS.

Required R packages are specified in the individual R scripts.

Package version information is provided in `package_versions.txt`.

No non-standard hardware is required.

## Installation

Install R and the R packages specified at the beginning of each script.

Typical installation time is approximately 10–20 minutes on a standard desktop computer.

## Usage

Run the corresponding R scripts after installing the required R packages.

Input file paths may need to be modified according to the local location of the input data.

## Demo

A small example dataset and R script are provided in the `demo` folder.

To run the demo, set the working directory to the root of this repository and run:

`source("demo/STEM121_demo.R")`

The demo performs Mann–Whitney U tests comparing TP and Combination groups at three distances, followed by Holm correction for multiple comparisons.

The expected output is `STEM121_demo_results.csv`, which will be generated in the `demo` folder.

The demo is expected to run within a few seconds on a standard desktop computer.
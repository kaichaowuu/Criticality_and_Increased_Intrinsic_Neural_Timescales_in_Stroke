# Criticality and Increased Intrinsic Neural Timescales in Stroke

This repository contains the analysis and modelling code accompanying our paper:

> Kaichao Wu, Beth Jelfs, Qiang Fang and Leonardo L. Gollo
> Criticality and increased intrinsic neural timescales in stroke
> npj Systems Biology and Applications 11.1 (2025): 103.

It includes Python, MATLAB, and Jupyter Notebook scripts used to compute intrinsic neural timescales (INT), their variants, dynamic measures, and the full modelling framework presented in the main text and supplementary materials.


If this repository is helpful, please cite us:

> @article{li2025machine,
> 
> title = {Criticality and increased intrinsic neural timescales in stroke},
> 
> author = {Kaichao Wu, Beth Jelfs, Qiang Fang and Leonardo L. Gollo}
> >
> journal = {npj Systems Biology and Applications},
> 
> volume = {12},
> 
> number = {4},
> 
> pages = {103},
> 
> year = {2025},
> 
> DOI = {https://doi.org/10.1038/s41540-025-00626-7},
> 
> publisher = {Nature Publishing Group UK London}
> }


## Repository Structure
├── Modelling

├── dynamic_measures.py

├── dynamic_measures.m

├── dynamic_measurement.py

├── Temporal_correlation.ipynb

├── stroke_int.ipynb

└── README.md

---

## Description of Files

### `dynamic_measures.py` & `dynamic_measurement.py`
Python scripts for computing intrinsic neural timescales (INT) and INT-related dynamic metrics, including HUST and other variants.  
Definitions and mathematical formulations of these variants can be found in our paper.

### `dynamic_measures.m`
MATLAB version of the intrinsic neural timescale computation pipeline.

### `stroke_int.ipynb`
Notebook for:
- Computing INT from fMRI data  
- Presenting stroke vs control comparisons  
- Reproducing empirical findings in the paper

### `Temporal_correlation.ipynb`
Extends the INT framework by analysing temporal correlation structures, complementing the results presented in the paper.

### `Modelling/`
Contains all modelling scripts used in the paper, including:
- Brain network modelling  
- Criticality analysis  
- Simulations based on structural connectivity  
- Results shown in both the main text and supplementary materials

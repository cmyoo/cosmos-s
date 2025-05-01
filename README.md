[![C++](https://custom-icon-badges.herokuapp.com/badge/C++-f34b7d.svg?logo=Cplusplus&logoColor=white)]()
[![Github All Releases](https://img.shields.io/github/downloads/atom/atom/total.svg?style=flat)]()
[![BSD-3-Clause](https://custom-icon-badges.herokuapp.com/badge/license-BSD%203%20Clause-8BB80A.svg?logo=law&logoColor=white)]()

# COSMOS-S
<!-- 
[![status](https://joss.theoj.org/papers/af52e7f1b7637bfa68818fde7c1a34de/status.svg)](https://joss.theoj.org/papers/af52e7f1b7637bfa68818fde7c1a34de)
[![DOI](https://zenodo.org/badge/118786602.svg)](https://zenodo.org/badge/latestdoi/118786602) -->

An open-source code for spherically symmetric numerical relativity simulation specialized for PBH formation. 
COSMOS is originally translated from SACRA code [Tetsuro Yamamoto, Masaru Shibata, Keisuke Taniguchi(arXiv:0806.4007)] into C++.
COSMOS-S is derived from COSMOS to be specialized for spherically symmetric situation.

Please visit https://sites.google.com/view/cosmoscode/ for a list of publications using COSMOS, and some others.

<!-- ## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRChombo/GRChombo/wiki), with the home page giving guidance on where to start. -->

## Documents
API documents are available in [https://cmyoo.github.io/cosmos-s/](https://cmyoo.github.io/cosmos-s/). 

Wiki pages are available in [https://github.com/cmyoo/cosmos-s/wiki](https://github.com/cmyoo/cosmos-s/wiki). 


## Description
Perfect fluid with linear equation of states & massless scalar field

Non-Cartesian scale-up coordinates

Fixed mesh refinement

OpenMP parallelization 

## Required system environments
Depending on the environment, one may have to edit the makefile appropriately. 

### Linux
the default makefile would work with g++ compiler.

### Mac
One may have to install an openMP package in addition to a C++ compiler. 
The std compiler option would have to be altered to an appropriate one (e.g. "-std=c++14")

<!-- ### Note
When performing test calculations on a low-spec PC, such as a laptop PC, the calculation speed (especially for the apparent horizon finder) may drop significantly unless the number of threads for openMP is kept small. For such test calculations, a user may have to limit the number of cores to one or at most a few with "export OMP_NUM_THREADS=x" (x: number of threads). -->

## Examples
Three samples are attached in the directories: **adiabatic_spherical_s**, **scalar_iso_s** and **typeIIB**. The detailed physical settings are given in the enclosed pdf files in the directories. 

## License
COSMOS-S is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Citation
If you use COSMOS-S as part of a paper, please cite the following three papers in which COSMOS and COSMOS-S are used:

```
@article{Yoo:2021fxs,
    author = "Yoo, Chul-Moon and Harada, Tomohiro and Hirano, Shin'ichi and Okawa, Hirotada and Sasaki, Misao",
    title = "{Primordial black hole formation from massless scalar isocurvature}",
    eprint = "2112.12335",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "YITP-21-161, RUP-21-23",
    doi = "10.1103/PhysRevD.105.103538",
    journal = "Phys. Rev. D",
    volume = "105",
    number = "10",
    pages = "103538",
    year = "2022"
}
@article{Okawa:2014nda,
    author = "Okawa, Hirotada and Witek, Helvi and Cardoso, Vitor",
    title = "{Black holes and fundamental fields in Numerical Relativity: initial data construction and evolution of bound states}",
    eprint = "1401.1548",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.89.104032",
    journal = "Phys. Rev. D",
    volume = "89",
    number = "10",
    pages = "104032",
    year = "2014"
}
@article{Yoo:2013yea,
    author = "Yoo, Chul-Moon and Okawa, Hirotada and Nakao, Ken-ichi",
    title = "{Black Hole Universe: Time Evolution}",
    eprint = "1306.1389",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    reportNumber = "OCU-PHYS:382, AP-GR:105, YITP-13-45",
    doi = "10.1103/PhysRevLett.111.161102",
    journal = "Phys. Rev. Lett.",
    volume = "111",
    pages = "161102",
    year = "2013"
}
```

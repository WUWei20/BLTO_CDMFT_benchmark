# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This repository is for benchmarking the cluster dynamical mean-field theory (CDMFT) study on the bilayer two-orbital Hubbard model of La3Ni2O7 [ YY Zheng and W Wú , arXiv:2312.03605 (2023)]
* Version: 0.1.0

### How do I get set up? ###

* You can run the open source DCA++ (https://github.com/CompFUSE/DCA/archive/refs/tags/paper.2019.new_code.zip)  to generate the benchmark data, with using C++ hearder files provided 
 in this repository that realizes the bilayer two-orbital Hubbard (BLTO) model in  DCA++. 
* Configuration: Please first download the DCA++ code (https://github.com/CompFUSE/DCA/archive/refs/tags/paper.2019.new_code.zip), then copy files in  "DCA++_additional_source_codes" folder
  to the root directory of DCA++ source code. Please keep the folder structure of "DCA++_additional_source_codes" unchanged!
* Dependencies: Just use the DCA++ dependencies.
* How to run tests: Use the DCA++ cmake script to compile. You can use following cmake variables (please see also DCA++ Docs): -DCMAKE_CXX_COMPILER=mpigxx -DDCA_WITH_TESTS_FAST=OFF -DTEST_RUNNER=mpirun -DCMAKE_BUILD_TYPE=Release -DDCA_BUILD_ANALYSIS=ON -DDCA_BUILD_DCA=ON -DDCA_POINT_GROUP=D4 -DDCA_WITH_MPI=ON -DDCA_LATTICE=BLTO -DDCA_RNG=std::mt19937_64 -DDCA_MODEL=tight-binding -DDCA_CLUSTER_SOLVER=CT-AUX

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

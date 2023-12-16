## Overview
  Demo code for VLDB paper :

  Catching More Associations by Referencing External Graphs

The code is developed based library GUNDAM https://github.com/MinovskySociety/GUNDAM

Main repo for graph rule discovery is being refactored and optimized, which will be released soon. https://github.com/MinovskySociety/GraphRules

## Install dependencies on Ubuntu

GCC version: 7.4.0 or above, support of c++17 standard required.

Install mpi:

```sh
sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
```
Install glog:

```sh
sudo apt-get install libgoogle-glog-dev
```

Install gflags:

```sh
sudo apt-get install libgflags-dev
```

Install yaml:

```sh
sudo apt-get install libyaml-cpp-dev
```


## Compile
```sh
mkdir build && cd ./build
cmake ../
make all -j
mv ./dataset/* ./build.
```
## Run
```sh 
srun -p $partition of machines -N $number of machines$ -n $number of processes -B ($number of CPUs):($number of cores per CPU):#(number of thread per core) ./gar_discover --yaml_file=./gar_discover.yaml
After running batch discovery, run C++ scripts in ./misc/* to build indices for incremental discovery.
srun -p $partition of machines -N $number of machines$ -n $number of processes -B ($number of CPUs):($number of cores per CPU):#(number of thread per core) ./inc_gar_discover --yaml_file=./inc_gar_discover.yaml
Note: (1) Please use at least two machines(for distributed setting) and at least two processes(for coordinator and workers).
      (2) The addresses of the dataset in yaml file of each process should be same.
      (3) Settings of batch and incremental discovery should be same.
```

## Details

```sh
We use libgrape-lite for multi-process parallelism and openmp for multi-thread parallelism;
Main lopp of batch discovery can be found in src/apps/rule_discover/gar_discover* .
Main loop of incremental discovery  can be found in src/apps/rule_discover/inc_gar_discover* .

Implementation of part important functions can be found in include/gar/* and src/apps/rule_discover/*
Some functions are shared by both batch and incremental discovery.
```

## Example
```sh
Examplary toy dataset is a fragment of dbpedia dataset, which describes scientific classification of around 10 thousands species.
Examplary rules (in ./dataset/gar*.csv) are rules discovered from the toy dataset.
```

```sh
For gar_discover.yaml,
(1) DataGraphPath is the path for data graph,
(2) KnowledgeGraphVFile and KnowledgeGraphEFile are addresses for knowledge graph,
(3) ERFile is the location of vertices mapping between data and knowledge graph,
(4) StoreMatch indicates whether the algorithm store indices for incremental discovery,
(5) MatchGARUsingSet indicates whether validates gars in group.
(6) PatternVertexLimit controls the maximum number of vertices of GARs.
 
For inc_gar_discover.yaml,
(1) SupportEstimation indicates whether the algorithm prunes GARs by estimating,
(2) EstimationThreshold is the threshold for estimation.
```

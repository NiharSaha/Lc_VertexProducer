# VertexCompositeAnalysis

Resonace decay reconstruction algorithms with ```VertexCompositeCandiate``` collection in cmssw. Compatible with 2023 PbPb datafomat. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

This branch was edited by Abigail Wesolek and currently will only process events of 

$D^{0} \to K+\pi$

But other branches can easily be added in the future.  You want to check the origianl branch for more information.

The code been edited ro un in 2 steps.
(1.) run_VCAnalyzer.py
(2.) run_VCProducer.py

Analyzer reads in a miniAOD file and outputs a skimmed EDM + VertexComositeCandidiate information (edm file format)
Producer reads in the previously produced edm file and outputs a root TTree.





## How to run

```bash 
#LXplus, bash, cmssw-el8 apptainer

cmsrel CMSSW_13_2_11

cd CMSSW_13_2_11/src
cmsenv
git cms-init

git clone git@github.com:a-wesole/VertexCompositeAnalysis.git
scram b -j8
cd VertexCompositeAnalysis/VertexCompositeProducer/test

cmsRun run_VCAnalyzer.py
cmsRun run_VCProducer.py

```


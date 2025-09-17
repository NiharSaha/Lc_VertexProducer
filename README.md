# VertexCompositeAnalysis

Resonace decay reconstruction algorithms with ```VertexCompositeCandiate``` collection in cmssw. Compatible with 2023 PbPb datafomat. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

This branch was edited by Nihar Ranjan Saha and currently will only process events of 

$\Lambda_{c}^{\pm} \to p^{\pm}K^{\mp}\pi^{\pm}$


The code has been edited to produce skimmedEDM and TTree simultaneously.


Analyzer reads in a miniAOD file and outputs a skimmed EDM + VertexComositeCandidiate information (edm file format)
Producer reads in the previously produced edm file and outputs a root TTree.





## How to run

#LXplus, bash, cmssw-el8 apptainer

cmsrel CMSSW_13_2_11
cd CMSSW_13_2_11/src
cmsenv

#clone the repo
git clone git@github.com:NiharSaha/Lc_VertexProducer.git

#run the setup.sh script to install the HeavyIonsAnalysis from CmsHI github, to be consistent with CMSSW version 
cd VertexCompositeAnalysis
./setup.sh

#compile
cd ..
scram b -j8

#cd and run code that produces edm and TTree files 
cd VertexCompositeAnalysis/VertexCompositeProducer/test
cmsRun Config_LcEdmTTree_Data.py  for Data
cmsRun Config_LcEdmTTree_MC.py    for MC



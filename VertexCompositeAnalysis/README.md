# VertexCompositeAnalysis

Resonace decay reconstruction algorithms with ```VertexCompositeCandiate``` collection in cmssw. Compatible with 2023 PbPb datafomat. The package is fully orthogonal to the ```HiForest``` framework to be combined with other objects.

This branch was edited by Nihar Ranjan Saha and currently will only process events of 

$\Lambda_{c}^{\pm} \to p^{\pm}K^{\mp}\pi^{\pm}$



Analyzer reads in a miniAOD file and outputs a skimmed EDM + VertexComositeCandidiate information (edm file format)
Producer reads in the previously produced edm file and outputs a root TTree.





## How to run




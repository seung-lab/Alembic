[![Build Status](https://travis-ci.org/seung-lab/Julimaps.svg?branch=master)](https://travis-ci.org/seung-lab/Julimaps)

# Alembic
ALignment of Electron Microscopy By Image Correlograms
A set of tools for elastic image registration in Julia.

## Alignment
Register coarsely-aligned sections to each other, and k-nearest neighboring sections, by blockmatching. Filter matches automatically and manually. Cover each section image in a triangle mesh, and globally deform all the meshes elastically to accommodate the correspondences. Render with a piecewise affine transform.

## Requirements  
```
julia v0.6
[cloud-volume](https://github.com/seung-lab/cloud-volume)
```

## Getting started  
Setup the Julia wrapper for CloudVolume by following the documentation [here](https://github.com/seung-lab/CloudVolume.jl).

Here's a simple startup script:
```
using Alembic
load_params(joinpath(homedir), ".julia/v0.6/Alembic/src/params/pinky40_test.json"); # specifies CloudVolume dirs and params
ms = make_stack(); # make meshset given the params specified in the JSON file
match!(ms); # blockmatch between the meshes in the meshset, per the params
elastic_solve!(ms); # relax the spring system
render(ms); # render the images and save to CloudVolume
```

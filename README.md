[![Build Status](https://travis-ci.org/seung-lab/Julimaps.svg?branch=master)](https://travis-ci.org/seung-lab/Julimaps)

# Alembic
ALignment of Electron Microscopy By Image Correlograms
A set of tools for elastic image registration in Julia.

# Terminology
* _Tile_: the base image unit from the microscope at the highest resolution
* _Section_: the image formed by combining all the tiles together
* _Overview_: a downsampled image of the entire substrate from which the tiles were imaged (a downsampled superset of the section)
* _Wafer_: a collection of sections, denoted by the number of substrates that can fit into the microscope at one time
* _Stack_: a collection of sections, that include all the wafers

# Pipeline
## Premontage
Register tile images to the overview image by cross correlating the entire tile (downsampled to the overview image resolution) across the entire overview image. Allow tile translations only. Review by inspecting overlay image of tiles on the overview. There is no intervention method, yet.
## Montage
Register tile images to one another by blockmatching in their overlapping portions, as calculated by the premontage translations. Cover each tile image in a triangle mesh, and deform the mesh elastically to accommodate the correspondences. Render with a piecewise affine transform. Review by inspecting the combined overlay plot of all tiles combined. Intervene by removing bad correspondences by identifying bad blockmatch images, or by clicking on correspondences that have non-consistent displacements relative to their neighbors.
## Prealignment
Register montaged sections to the previous prealigned section by blockmatching. Render with a regularized affine transform (part affine, part rigid). Review by looking at overlays of the two sections. Intervene by removing bad correspondences by identifying bad blockmatch images, or by clicking on correspondences that have non-consistent displacements relative to their neighbors. 
## Alignment
Register prealigned sections to each other, and N+k neighboring sections, by blockmatching. Cover each section image in a triangle mesh, and globally deform all the meshes elastically to accommodate the correspondences. Render with a piecewise affine transform. Review by inspecting the combined overlay plot of all tiles combined. Intervene by removing bad correspondences by identifying bad blockmatch images, or by clicking on correspondences that have non-consistent displacements relative to their neighbors.


# Getting started  
Alembic.jl requires Julia v0.6, and depends on a chunked-volume storage system that can be interfaced with through [CloudVolume](https://github.com/seung-lab/cloud-volume).  

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

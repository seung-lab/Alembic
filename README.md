# Julimaps
A set of tools for elastic image registration in Julia.

# To do
* Find triangle that a point belongs to
* Finding max of cross correlation with subpixel accuracy
* Finding eignevalue ratio for principal curvature
* Try VoronoiDelaunay for mesh generation
* Check how S&co generate their spring mesh
* Create imfuse equivalent for rendering multiple images as one
  * this could be helpful for blending https://github.com/timholy/Images.jl/blob/master/doc/overlays.md
* Finish imwarp
* Optimize incidences2triangles
* Optimize pa_warp2 (rectify number of operations with time spent in method)
* Visualize block matching locations
* Test better interpolation method in piecewise affine warping
* Detect overlapping tile pairs based on affine transforms
* Include unit tests in MeshSolve
* Include unit tests in Mesh
* Include unit tests in piecewiseaffine_warp
* Create downsampling method
* Test imfuse on stitching more than two tiles

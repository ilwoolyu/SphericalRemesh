# SphericalRemesh: Sphere Mesh and Data Manipulation

## Description
This tool supports sphere mesh manipulation including rigid/non-rigid deformation as well as re-tessellation of the input sphere. This tool has been developed to support [HSD](https://github.com/ilwoolyu/HSD) that offers decomposable spherical deformation. In [HSD](https://github.com/ilwoolyu/HSD), the spherical registration can be achieved by fitting spherical harmonics coefficients. This tool takes spherical harmonics coefficients to deform the input (moving) sphere mesh and help adjust the number of the basis functions to control a level of deformation, which is useful for <i>**spherical data augmentation**</i>. The output spherical data can also be re-tessellated regarding a reference sphere mesh, by which the output data have the same structure as that of the reference sphere mesh.

## Installation
You can download and compile the source code using <a href="https://cmake.org/">CMake</a>. You can also download the [singularity image](https://vanderbilt.box.com/v/cmorph-release) or pull <a href="https://hub.docker.com/r/ilwoolyu/cmorph/">docker image</a>:
```
$ docker pull ilwoolyu/cmorph:<version>
```
Click [here](https://hub.docker.com/r/ilwoolyu/cmorph) to check how to use the tool.

## Requirements for build
* [MeshLib](https://github.com/ilwoolyu/MeshLib)

## Usage
* Inputs
  * Sphere mesh (`vtk`)
  * Spherical data (`txt`)
  * Reference sphere mesh (`vtk`)
  * (optional) Spherical harmonics coefficients (`txt`)
* Outputs
  * Spherical data on the reference sphere mesh (`txt`)
  * Deformed spherical data (`txt`) and/or deformed sphere mesh (`vtk`) if spherical harmonics are provided.

## Re-tessellation of Spherical Data
To deform data after HSD registration (see how to register spheres using HSD in a [group-wise](https://github.com/ilwoolyu/HSD#usage) or [pair-wise](https://github.com/ilwoolyu/HSD#pairwise-registration) manner):
```
$ SphericalRemesh \
-s <input_sphere> \
-r <reference_sphere> \
-c <coefficients_from_HSD> \
-p <a list of input_spherical_data> \
--outputProperty <output_prefix_for_spherical_data>
```
For example, the command below generates deformed spherical data tessellated by `ico7.vtk`.
```
$ SphericalRemesh \
-s lh.sphere.vtk \
-r ico7.vtk \
-c lh.coeff.txt \
-p lh.curv.txt lh.sulc.txt \
--outputProperty lh.reg
```
The following outputs will be generated:
```
$ ls lh.reg.*
lh.reg.curv.txt lh.reg.sulc.txt
```
If you have a registered sphere mesh, the following command is equivalent to the above one:
```
$ SphericalRemesh \
-s lh.sphere.reg.vtk \
-r ico7.vtk \
-p lh.curv.txt lh.sulc.txt \
--outputProperty lh.reg
```
Similarly, if you need only re-tessellated spherical data, ignore `-c`. This is useful to create a consistent mesh structure across different sphere mesh files, which, of course, does not require spherical harmonics coefficients.

If you only need a sphere mesh deformed by HSD, use the following command:
```
$ SphericalRemesh \
-s <input_sphere> \
-c <coefficients_from_HSD> \
--deform <output_sphere>
```

## Fast Spherical Data Augmentation
For data augmentation [[2](#ref2),[3](#ref3)], you can set a range of spherical harmonics using `--deg0` and `--deg` flags. By adding the two flags to the same example above, you have
```
$ SphericalRemesh \
-s lh.sphere.vtk \
-r ico7.vtk \
-c lh.coeff.txt \
-p lh.curv.txt lh.sulc.txt \
--outputProperty lh.reg.aug \
--deg0 0
--deg 3
```
This will augment spherical data.
```
$ ls lh.reg.aug*
lh.reg.aug0.curv.txt lh.reg.aug0.sulc.txt lh.reg.aug1.curv.txt lh.reg.aug1.sulc.txt
lh.reg.aug2.curv.txt lh.reg.aug2.sulc.txt lh.reg.aug3.curv.txt lh.reg.aug3.sulc.txt
```
>**Note1** `aug#` is a level of deformation. `aug0` captures rigid rotation and the deformation becomes more non-rigid and closer to the target as `#` increases. `#` supports up to the maximum degree of spherical harmonics. See `-d` flag in [HSD](https://github.com/ilwoolyu/HSD#usage).

>**Note2** For the categorical data, add `--nneighbor` to avoid barycentric interpolation of the data. The tool uses barycentric interpolation for spherical data re-tesselation by default, otherwise.

## Re-tessellation for Generic Mesh
You can re-tessellate a generic mesh if its associated sphere mesh is available. Use the following command:
```
$ SphericalRemesh \
-i <input_surface>
-s <input_sphere> \
-r <reference_sphere> \
-o <output_surface> \
-c <coefficients_from_HSD> # optional
```
The output mesh will have the same structure as that of the reference sphere.

## References
Please cite the following papers if you find it useful.
* <a id="ref1"></a>Lyu, I., Kang, H., Woodward, N., Styner, M., Landman, B., <a href="https://doi.org/10.1016/j.media.2019.06.013">Hierarchical Spherical Deformation for Cortical Surface Registration</a>, <i>Medical Image Analysis</i>, 57, 72-88, 2019</li>
* <a id="ref4"></a>Lyu I., Bao, S., Hao, L., Yao, J., Miller, J., Voorhies, W., Taylor, W., Bunge, S., Weiner, K., Landman, B., <a href="https://doi.org/10.1016/j.neuroimage.2021.117758">Labeling Lateral Prefrontal Sulci using Spherical Data Augmentation and Context-aware Training</a>, <i>NeuroImage</i>, 229, 117758, 2021

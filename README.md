# Stitch Meshing

[Kui Wu](http://www.cs.utah.edu/~kwu/), 
[Xifeng Gao](https://gaoxifeng.github.io/), 
[Zachary Ferguson](http://zfergus.me/), 
[Daniele Panozzo](https://cs.nyu.edu/~panozzo), 
[Cem Yuksel](http://www.cemyuksel.com/)

*Siggraph 2018*

## Abstract

We introduce the first fully automatic pipeline to convert arbitrary 3D shapes into knit models. Our pipeline is based on a global parametrization remeshing pipeline to produce an isotropic quad-dominant mesh aligned with a 2-RoSy field. The knitting directions over the surface are determined using a set of custom topological operations and a two-step global optimization that minimizes the number of irregularities. The resulting mesh is converted into a valid stitch mesh that represents the knit model. The yarn curves are generated from the stitch mesh and the final yarn geometry is computed using a yarn-level relaxation process. Thus, we produce topologically valid models that can be used with a yarn-level simulation. We validate our algorithm by automatically generating knit models from complex 3D shapes and processing over a hundred models with various shapes without any user input or parameter tuning. We also demonstrate applications of our approach for custom knit model generation using fabrication via 3D printing.

 ## Project Website
 
 http://www.cs.utah.edu/~kwu/stitchmodeling#stitchmeshing

## Installation
- Clone the repository into your local machine:

```bash
git clone https://github.com/kuiwuchn/stitchMeshing --recursive
```
- Compile the code using cmake 

You need to install [Gurobi](http://www.gurobi.com/) before compiling the code.

Set include directory and lib directory accordingly for gurobi in CMakeLists.txt

```
cd stitchMeshing
mkdir build
cd build
cmake ..
make
```

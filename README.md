# FiniteMesh.jl

[![version](https://juliahub.com/docs/FiniteMesh/version.svg)](https://juliahub.com/ui/Packages/FiniteMesh/zdt25)
[![CI](https://github.com/vavrines/MeshArt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/vavrines/MeshArt.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/vavrines/FiniteMesh.jl/branch/main/graph/badge.svg?token=DpQ1YKKlBL)](https://codecov.io/gh/vavrines/FiniteMesh.jl)
[![deps](https://juliahub.com/docs/FiniteMesh/deps.svg)](https://juliahub.com/ui/Packages/FiniteMesh/zdt25?t=2)

This package provides lightweight methods for mesh I/O. 
The input interface is based on the Python project [meshio](https://github.com/nschloe/meshio).
The supported formats include:
> [Abaqus](http://abaqus.software.polimi.it/v6.14/index.html) (`.inp`),
 ANSYS msh (`.msh`),
 [AVS-UCD](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) (`.avs`),
 [CGNS](https://cgns.github.io/) (`.cgns`),
 [DOLFIN XML](https://manpages.ubuntu.com/manpages/disco/man1/dolfin-convert.1.html) (`.xml`),
 [Exodus](https://nschloe.github.io/meshio/exodus.pdf) (`.e`, `.exo`),
 [FLAC3D](https://www.itascacg.com/software/flac3d) (`.f3grid`),
 [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) (`.h5m`),
 [Kratos/MDPA](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) (`.mdpa`),
 [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) (`.mesh`, `.meshb`),
 [MED/Salome](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) (`.med`),
 [Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) (bulk data, `.bdf`, `.fem`, `.nas`),
 [Neuroglancer precomputed format](https://github.com/google/neuroglancer/tree/master/src/neuroglancer/datasource/precomputed#mesh-representation-of-segmented-object-surfaces),
 [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) (format versions 2.2, 4.0, and 4.1, `.msh`),
 [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) (`.obj`),
 [OFF](https://segeval.cs.princeton.edu/public/off_format.html) (`.off`),
 [PERMAS](https://www.intes.de) (`.post`, `.post.gz`, `.dato`, `.dato.gz`),
 [PLY](https://en.wikipedia.org/wiki/PLY_(file_format)) (`.ply`),
 [STL](https://en.wikipedia.org/wiki/STL_(file_format)) (`.stl`),
 [Tecplot](http://paulbourke.net/dataformats/tp/)(`.dat`),
 [TetGen](https://wias-berlin.de/software/tetgen/fformats.html)(`.node/.ele`),
 [SVG](https://www.w3.org/TR/SVG/) (2D output only) (`.svg`),
 [SU2](https://su2code.github.io/docs_v7/Mesh-File) (`.su2`),
 [UGRID](http://www.simcenter.msstate.edu/software/downloads/doc/ug_io/3d_grid_file_type_ugrid.html) (`.ugrid`),
 [VTK](https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf) (`.vtk`),
 [VTU](https://www.vtk.org/Wiki/VTK_XML_Formats) (`.vtu`),
 [WKT](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) ([TIN](https://en.wikipedia.org/wiki/Triangulated_irregular_network)) (`.wkt`),
 [XDMF](https://www.xdmf.org/index.php/XDMF_Model_and_Format) (`.xdmf`, `.xmf`).

## Installation

FiniteMesh.jl is a registered package in the official Julia package registry. 
We recommend installing it with the built-in Julia package manager.
From the Julia REPL, you can get in the package manager (by pressing ]) and add the package.
This will automatically install the package and all its dependencies.

```julia
julia> ]
(v1.6) pkg> add FiniteMesh
```

## Usage

To read a mesh, simply do
```julia
using FiniteMesh
cells, points = read_mesh("path-of-mesh-file")
```
The resulted `points` are the coordinates of nodes, and `cells` provides the affiliation information of these points to cell IDs.

The connectivity calculator is provided with native Julia.
The following information can be inferred from `cells` and `points`: 
- Node IDs of a face
- Adjacent cell IDs of a face
- Neighbor cell IDs of a cell
- Cell type (inner / boundary)
- Cell volume
- Face area
- Cell center location
- Face center location
- Unit normal vector of a face

For convenience, a `Mesh` struct is defined to contain all the above information as its fields.
It can be constructed simply via
```julia
mesh = Mesh("path-of-mesh-file")
```
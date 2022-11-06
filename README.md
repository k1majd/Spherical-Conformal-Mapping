# Spherical Conformal Mapping

This package performs Spherical Conformal Mapping on a given triangular mesh saved in .off mesh format.
For the construction of triangular mesh from a given .off file, the **halfedge_mesh** package in 
https://github.com/carlosrojas/halfedge_mesh.git is used. 

# Usage
To perform the mapping: 
1. Clone this repository.
```
git clone https://github.com/k1majd/Spherical-Conformal-Mapping.git
```
2. Store your .off file in `halfedge_mesh/tests/data/`
3. Run `main.py -m [name of mesh file in string form without .off]` in `halfedge_mesh/`. For example:
```
main.py -m "brain" 
```

# Manual

## Config file
The config file has the following options:
* Material Constants:
  * damping constant
  * magnetic moment
  * uniaxial anisotropy
  * lattice parameter  
* Lattice Vectors
* Spin sites within unit cell (as a fraction of the lattice vectors)
* The initial Magnetisation of the sites
* Sites as integer positions (For spinwave calculations)
* External Applied field:
  * Field type. Current options include _uniform, split, sine_.
  * Field size (for _split_ and _uniform_)
* Time settings:
  * size of timestep in seconds
  * Number of timesteps
  * Relaxation time (for calculating equilibrium spin configuration)
* System Geometry:
  * System size:
    * number of units cells in x direction.
    * number of units cells in y direction.
    * number of units cells in z direction.
  * Number of sites within the unit cell
  * Number of sublattices within the system.

## Compilation

## Cuda

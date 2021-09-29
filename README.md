
2D TDSE Solver
===============================
# About

A two-dimensional time-dependent Schrödinger equation visualiser that illustrates the dynamics of a gaussian wavepacket scattered off of different potential fields.

Numerical solutions to the TDSE can be obtained by using an ADI or Pseudospectral FFT based method.


# Usage

Using Ctrl+a and Ctrl+Enter in VSCode from "main.jl" is the recommended way to run the code. Does NOT work directly from terminal ( This might just be my system though ). 



User interaction is enabled for specifying the following,

### Domain parameters
```
N     : Number of grid points along each axis ( N = Nx = Ny )
Nt    : Number of time steps
t_max : Simulated runtime
x, x  : Boundary points of domain along x-axis
y, y  : Boundary points of domain along y-axis
```

### Wave parameters
```
σx, σy   : Initial spatial spread 
kx, ky   : Initial wavenumber
x_0, y_0 : Coordinates wave packet is centered at 
```

### Potential Field 
```
Potential well

Barrier
    - pos   : position double slit barrier is centered at
    - gap   : slit size
    - width : width of barrier

Cylinder 
    - x, y  : position potential is centered at
    - r     : radius of circle/cylinder
    - |U|   : magnitude of potential field
```
# Examples
![diffraction_resize](https://user-images.githubusercontent.com/81137805/135265939-afa2d024-cc1e-4ed9-8fe6-c36559cac022.gif)
![cylinder_resize0](https://user-images.githubusercontent.com/81137805/135265989-486c3b8c-87eb-4c75-91d3-a571d49327b9.gif)
![cylinder_resize](https://user-images.githubusercontent.com/81137805/135271895-93a4f0ad-ebc3-4964-94f7-9272d69bb02b.gif)
# Version Info
```
Julia   v1.6.2
GLMakie v0.4.4
FFTW    v1.4.3
```

# References for ADI solver

```
Galbraith I, Ching Y, Abraham E. Two‐dimensional time‐dependent quantum‐mechanical scattering event. American Journal of Physics. 1984;52(1):60-68. 
```






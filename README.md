# VortexLattice.jl

VortexLattice.jl is a Julia package for the Vortex Lattice Method I've created for my hobbies and academic projects.
It is meant to make simple the conceptual design of aircraft configurations of interest, as in the example below.

It is currently capable of obtaining:

* Hinge moments;
* Sectional and surface coefficients;
* Full configuration coefficients;
* Geometry and spatial distribution plots (using matplotlib as a back-end).

## Defining Viscous Corrections for Airfoils

Struct to contain information about an airfoil:

* `camberline_coords`: matrix `(shape (npts, 2))` containing airfoil camber line coordinates
* `CL_correct`: function as `C_L=C_L_inv+CL_correct(CL)`
* `CD_correct`: function as `C_D=C_D_inv+CD_correct(CL)`
```
mutable struct Airfoil
    camberline_coords::Matrix{Float64}
    CL_correct::Function
    CD_correct::Function
end
```

Constructor for inviscid flat plate airfoil (default for wing sections):
```
inviscid_flatplate()=Airfoil(reshape([0.0, 1.0, 0.0, 0.0], 2, 2), CL->0.0, CL->0.0)
```

## Defining Lifting Surfaces

Constructor for wing section (uses a flat plate airfoil as default):

* `LE`: position of the leading edge
* `chord`: section chord
* `incidence`: section incidence. Defaults to zero
* `afl`: instance of struct Airfoil to use for viscous corrections. Defaults to flat plate
* `controls`: array of tuples indicating present control surfaces. Tuples should indicate, respectively, the control's name (Symbol), the percentage of the chord in
which it's initiated and its gain in the present section
```
Section(LE::Vector{Float64}, chord::Float64; 
    incidence::Float64=0.0, afl=inviscid_flatplate(), 
    controls=[])=Section(LE, chord, incidence, afl, controls)
```

Constructor for a lifting surface based on a list of sections (from left to right):

For arguments `cdiscs` and `bdiscs`, a usage example is:

```
sects=[
    (LE=[1.0, -1.0, 0.0], chord=0.5),
    (LE=[0.0, 0.0, 0.0], chord=1.0),
    (LE=[1.0, -1.0, 0.0], chord=0.5)
]
bdiscs=[
    [0.0, 0.5, 1.0],
    [0.0, 1.0]
]
cdiscs=[
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0],
    [0.0, 0.5, 1.0]
]
```

Which would result in the panel corners:
```
[
    [1.0, -1.0, 0.0] [0.5, -0.5, 0.0] [0.0, 0.0, 0.0] [1.0, 1.0, 0.0]
    [1.25, -1.0, 0.0] [0.875, -0.5, 0.0] [0.5, 0.0, 0.0] [1.25, 1.0, 0.0]
    [1.5, -1.0, 0.0] [1.25, -0.5, 0.0] [1.0, 0.0, 0.0] [1.5, 1.0, 0.0]
]
```

* `sects`: vector of sections to base the surface on
* `cdiscs`: vector of vectors, each describing a mapping for the chordwise discretization in a wing section (ranging from zero to 1, and spaced as desired for panel
edges)
* `bdiscs`: vector of vectors, each describing a mapping for the spanwise discretization in a wing quadrant
```
function Surface(sects::Vector{Section}, cdiscs::Vector{Vector{Float64}}, 
    bdiscs::Vector{Vector{Float64}})
```

## Simulating an Aircraft

Constructor for Aircraft struct based on an array of surfaces:

* `surfaces`: vector of surfaces
* `Sref`: reference surface
* `cref`: reference chord
* `bref`: reference span

```
function Aircraft(surfaces::Vector{Surface}; 
    Sref::Float64=1.0, cref::Float64=1.0, bref::Float64=1.0)
```

Function to obtain all pertinent coefficients for an aircraft:

* `acft`: `Aircraft` struct instance
* `alpha`: angle of attack
* `beta`: sideslip angle
* `p`: x axis (stability and control coordinate system) angular velocity
* `q`: y axis (stability and control coordinate system) angular velocity
* `r`: z axis (stability and control coordinate system) angular velocity
* `control_deflections`: dictionary relating control symbols and control deflections
* `CG`: vector with coordinates for momentum calculations
```
function get_data(acft::Aircraft; alpha::Float64=0.0, beta::Float64=0.0,
    p::Float64=0.0, q::Float64=0.0, r::Float64=0.0, control_deflections::Dict{Symbol, Float64}, 
    CG::Vector{Float64}=zeros(Float64, 3))
```

`get_data` should return a dictionary with keys for sectional forces, surface forces, total forces and hinge moments. An example of the output data structure for the output dictionary can be seen by running `examples/tapered_swept_wing.jl`.

To obtain the distribution of positions and chords around which sectional data is produced, use vectors `Surface.quarter_chords` and `Surface.chords`.

## Further reference

For more details, refer to the docstrings with Julia's help menu, or locally generate the automatic documentation with:

```
$ julia
include("docs/make.jl")
exit()
$ firefox docs/build/index.html
```

## Installation

To install VortexLattice.jl, use:

```
add https://github.com/pedrosecchi67/VortexLattice
```

... in the Julia package manager.

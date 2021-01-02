# VortexLattice.jl

VortexLattice.jl is a Julia package for the Vortex Lattice Method I've created for my hobbies and academic projects.
It is meant to make simple the conceptual design of aircraft configurations of interest, as in the example below.

It is currently capable of obtaining:

* Hinge moments;
* Sectional and surface coefficients;
* Full configuration coefficients;
* Geometry and spatial distribution plots (using matplotlib as a back-end).

```
using VortexLattice

# geometry

b=1.0
croot=0.25
taper=0.5
twist=-10.0
theta_root=5.0
sweep=30.0
dihedron=5.0

# flight conditions

alpha=5.0
delta_aileron=5.0

# reference dimensions

Sref=b*croot*(1.0+taper)/2
cref=(2.0*(taper^2+taper+1.0))*croot/(3.0*(1.0+taper))
bref=b

# discretization settings

cos_dist=collect(LinRange(-pi/2, pi/2, 20))
cos_dist=[sin(e)/2+0.5 for e in cos_dist]

# airfoils (see API reference for further detail)

fplate=inviscid_flatplate()

# wing sections

sleft=Section(
    [b*tand(sweep)/2, -b/2, b*tand(dihedron)/2], 
    croot*taper; 
    incidence=theta_root+twist, 
    controls=[(:aileron, 0.2, 1.0)]
)
scenter=Section(
    [0.0, 0.0, 0.0], 
    croot; incidence=theta_root, 
    controls=[(:aileron, 0.2, 1.0), 
    (:aileron, 0.2, -1.0)]
)
sright=Section(
    [b*tand(sweep)/2, b/2, b*tand(dihedron)/2], 
    croot*taper; 
    incidence=theta_root+twist, 
    controls=[(:aileron, 0.2, -1.0)]
)

cdiscs=[cos_dist, cos_dist, cos_dist]
bdiscs=[cos_dist, cos_dist]

# surfaces

surf=Surface([sleft, scenter, sright], cdiscs, bdiscs)

# complete aircraft and solution

acft=Aircraft([surf]; Sref=Sref, bref=bref, cref=cref)

# coefficient calculations

dat=get_data(acft; alpha=alpha, control_deflections=Dict(:aileron=>delta_aileron))

# plotting (based on matplotlib wrappers)

plot_aircraft(acft)
```

## Installation

To install VortexLattice.jl, use:

```
add https://github.com/pedrosecchi67/VortexLattice
```

... in the Julia package manager.

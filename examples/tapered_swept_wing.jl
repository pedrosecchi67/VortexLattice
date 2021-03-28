using VortexLattice

b=1.0
croot=0.25
taper=0.5
twist=-10.0
theta_root=5.0
sweep=30.0
dihedron=5.0

alpha=5.0
delta_aileron=5.0

Sref=b*croot*(1.0+taper)/2
cref=(2.0*(taper^2+taper+1.0))*croot/(3.0*(1.0+taper))
bref=b

cos_dist=cossenoidal_distribution(20)

afl = Airfoil(
    "gehesky01.dat",
    0.9191,
    [(0.5290, 0.0594), (0.9766, 0.0151), (2.1697, 0.0378)]
)

sleft=Section([b*tand(sweep)/2, -b/2, b*tand(dihedron)/2], croot*taper; incidence=theta_root+twist, controls=[(:aileron, 0.2, 1.0)], afl = afl)
scenter=Section([0.0, 0.0, 0.0], croot; incidence=theta_root, controls=[(:aileron, 0.2, 1.0), (:aileron, 0.2, -1.0)], afl = afl)
sright=Section([b*tand(sweep)/2, b/2, b*tand(dihedron)/2], croot*taper; incidence=theta_root+twist, controls=[(:aileron, 0.2, -1.0)], afl = afl)

cdiscs=[cos_dist, cos_dist, cos_dist]
bdiscs=[cos_dist, cos_dist]

surf=Surface([sleft, scenter, sright], cdiscs, bdiscs)

acft=Aircraft([surf]; Sref=Sref, bref=bref, cref=cref)

dat=get_data(acft; alpha=alpha, control_deflections=Dict(:aileron=>delta_aileron))

# plot_aircraft(acft)

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

cos_dist=collect(LinRange(-pi/2, pi/2, 20))
cos_dist=[sin(e)/2+0.5 for e in cos_dist]

fplate=inviscid_flatplate()

sleft=Section([b*tand(sweep)/2, -b/2, b*tand(dihedron)/2], croot*taper; incidence=theta_root+twist, controls=[(:aileron, 0.2, 1.0)])
scenter=Section([0.0, 0.0, 0.0], croot; incidence=theta_root, controls=[(:aileron, 0.2, 1.0), (:aileron, 0.2, -1.0)])
sright=Section([b*tand(sweep)/2, b/2, b*tand(dihedron)/2], croot*taper; incidence=theta_root+twist, controls=[(:aileron, 0.2, -1.0)])

cdiscs=[cos_dist, cos_dist, cos_dist]
bdiscs=[cos_dist, cos_dist]

surf=Surface([sleft, scenter, sright], cdiscs, bdiscs)

@time begin
  acft=Aircraft([surf]; Sref=Sref, bref=bref, cref=cref)
end

@time begin
  dat=get_data(acft; alpha=alpha, control_deflections=Dict(:aileron=>delta_aileron))
end

plot_aircraft(acft)

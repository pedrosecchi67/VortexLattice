using VortexLattice
using Plots

b = 3.0
croot = 0.6363
taper = 0.449
twist = -2.0
theta_root = 0.0
sweep = -19.0
dihedron = 2.5

delta_aileron = 0.0
h = 0.2 # altura da aeronave 

Sref = b * croot * (1.0 + taper) / 2
cref = (2.0*(taper^2 + taper + 1.0))*croot/(3.0*(1.0+taper))
bref = b

cos_dist = cossenoidal_distribution(20)

afl = Airfoil(
    "alexky.dat",
    1.03,
    [(0.1370, 0.01), (0.2303, 0.00833), (1.2768, 0.01558)]
)

sleft = Section([b*tand(sweep)/2, -b/2, b*tand(dihedron)/2], croot*taper;
                incidence=theta_root+twist, controls=[(:aileron, 0.2, 1.0)], afl=afl)
scenter = Section([0.0, 0.0, 0.0], croot;
                  incidence=theta_root, controls=[(:aileron, 0.2, 1.0), (:aileron, 0.2, -1.0)], afl=afl)
sright = Section([b*tand(sweep)/2, b/2, b*tand(dihedron)/2], croot*taper;
                 incidence=theta_root+twist, controls=[(:aileron, 0.2, -1.0)], afl=afl)

cdiscs = [cos_dist, cos_dist, cos_dist]
bdiscs = [cos_dist, cos_dist]

surf = VortexLattice.Surface([sleft, scenter, sright], cdiscs, bdiscs)
acft = Aircraft([surf]; Sref=Sref, bref=bref, cref=cref)
acft_ge = Aircraft([surf]; Sref=Sref, bref=bref, cref=cref, ground_effect=true, ground_height=h)

alphas = range(-10, 15, length=26)

# Teste sem efeito solo (GE desativado)
dat = [get_data(acft; alpha=alpha, control_deflections=Dict(:aileron => delta_aileron)) for alpha in alphas]

# Teste com efeito solo ativado (GE habilitado), altura 
dat_ge = [get_data(acft_ge; alpha=alpha, control_deflections=Dict(:aileron => delta_aileron)) for alpha in alphas]

# Plotando o gráfico de CL vs alpha, CD vs alpha e CL x CD com e sem o efeito solo
p1 = plot(alphas, [dat[i][:total][:CL] for i in eachindex(alphas)], xlabel="α [deg]", ylabel="CL", title="CL x α", label="Sem GE")
plot!(alphas, [dat_ge[i][:total][:CL] for i in eachindex(alphas)], label="No GE")
xlims!(-10, 15)
ylims!(-1.0, 2.5)

p2 = plot(alphas, [dat[i][:total][:CD] for i in eachindex(alphas)], xlabel="α [deg]", ylabel="CD", title="CD x α", label="Sem GE")
plot!(alphas, [dat_ge[i][:total][:CD] for i in eachindex(alphas)], label="No GE")
xlims!(-10, 15)
ylims!(0, 0.2)

p3 = plot([dat[i][:total][:CD] for i in eachindex(alphas)], [dat[i][:total][:CL] for i in eachindex(alphas)], xlabel="CD", ylabel="CL", title="CL vs CD", label="Sem GE")
plot!([dat[i][:total][:CD] for i in eachindex(alphas)], [dat_ge[i][:total][:CL] for i in eachindex(alphas)], label="No GE")

layout = @layout [a; b c]
graph = plot(p3, p1, p2, layout=layout)
display(graph)

# Plotando aeronave
plot_aircraft(acft)

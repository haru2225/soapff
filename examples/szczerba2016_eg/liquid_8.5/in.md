# liquid ethylene glycol, OPLS-AA (soapff oplsaa_eg), 293 K, 1 bar
include ff/system.settings
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 1000
thermo_style custom step temp press pe density
min_style cg
minimize 1.0e-4 1.0e-6 2000 20000
reset_timestep 0
velocity all create 293 5 dist gaussian mom yes rot yes
timestep 1.0
fix 1 all nvt temp 293 293 100
run 100000
unfix 1
reset_timestep 0
fix 2 all npt temp 293 293 100 iso 1 1 1000
variable rho equal density
fix 3 all ave/time 100 10 1000 v_rho file density.dat
run 400000
write_data final.data

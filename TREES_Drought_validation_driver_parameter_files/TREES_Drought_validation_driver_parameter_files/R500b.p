2180.0	altitude, m@Larm
41.367 	latitude
-105.238	longitude
0.5	z_ref, m 
2.24105406 0.1 20.0	lai, single sided
1.0 0.1 1.00  canopy_cover
0.1	canopy_height, m
0.00	fcloud, may model
1.0	l_angle, spherical, may sample
0.97	canopy emissivity
0.5	fPAR_beam
0.5	fPAR_diff
0.8	alpha_PAR
0.2	alpha_NIR
1.0 0.5 1.00	omega
2.0	p_crown
0.5	d_factor, C&Nfig5.5
0.1	zm_factor, C&Nfig5.5
0.2	zh_factor, 7.19
.0292 0.001 0.010      Rd_mult, Rd=Rd_mult*Vcmax
0.16 0.014 3.0    regression slope of phi2_300 and leaf water potential (Jmax_mult originally)
0.55 0.7 0.99     saturated phi2 at 300 mol PAR (originally thetaJ)
0.77 0.01 0.5   alpha psII quantum efficiency as PAR approaches zero (originally phij_sun)
0.04 0.01 0.5   kappa psII quantum efficiency as PAR approaches infinity (originally phiJ_shd)
0.0011 0.0012 0.0900    Nleaf
1.0 0.78513 1.0 N_fixed_proportion
0.2 0.01 0.2  Nrubisco, proportion
38.67764 23.9 62.9  Kc25, (Pa) MM const carboxylase, 25 deg C
2.1    q10Kc, (DIM) Q_10 for kc (default 2.1)
26123.26 17600 51600   Ko25, (Pa) MM const oxygenase, 25 deg C
1.2    q1 Ko, (DIM) Q_10 for ko
3.6    act25, (umol/mgRubisco/min) Rubisco activity
2.4 2.0 3.0    q10act, (DIM) Q_10 for Rubisco activity (default was 2.4)
0.24 0.86 2.65 	Gsref
0.3 0.47  	m (proportion of Gsref)
-0.05   Md
-0.34 -0.105 -0.18 -1.6 midday_at_sat_kl
6.5 7.7 0.5 3.0 e_at_saturated_kl
4.0 rhizosphere_width_(mm)
0.001 E_inc(Sperry98_pg351_last_para.)
4 soilshells
0.18 GMP_(mm)_geometric_mean_particle_diameter
8.5 GSD_geometric_standard_deviation_of_particle_size
0.8 BD_(MG/m3)_soil_bulk_density
0.39  0.44 0.5 porosity
0.43 silt_fraction
0.23 clay_fraction
1.0 frac_absorbing_length
0.01 0.1 10.0 Capacitance_(mol/Mpa*m2)_on_leaf_area_basis
1.0 axK:latKr_shoot_modules
1.0 axkr:latKr_root_modules
50.0 %total_R_in_root_system
18.0 saturated_kl_for_whole_plant (recalculated in the model)
-0.1 -0.1 0.05  pd_at_sat_kl
2.35 1.57 2 0.9 1.5 ax_Shoot-b_value_(weibull)
2.0 0.9 1.5 ax_Shoot-c_value_(weibull)
2.35 1.57 0.9 1.5 lat_Shoot-b_value_(weibull)
2.0 0.9 1.5 lat_Shoot-c_value_(weibull)
2.35 1.57 1.5 ax_Root-b_value_(weibull)
2.0 1.5 ax_Root-c_value_(weibull)
2.35 1.57 1.383 1.5 lat_Root-b_value_(weibull)
2.0 1.48 1.5 lat_Root-c_value_(weibull)
60.0 60.0 initial_conductivity(root)
0.01 decrement(root)- default 0.001
80.0 80.0 initial_conductivity(shoot)
0.1 decrement(shoot)
0.18 0.05 0.48 theta_opt
30.0 25.0 45.0 optimal_soil_T
1.0   growth_resp_proportion
0.0011 resp_coef_root, kg kg-1 day-1 deg 
0.0002 resp_coefficient_stem, kg kg-1 day-1 deg
0.0004 resp_coefficient_leaf, kg kg-1 day-1 deg
0.05 0.085 resp_coefficient (Q10), degC-1
72.26 71.22 73.30 EaSx, kjmol-1
0.000000995 0.000000877 0.00000111 kMsx, gCcm-3soil
0.000000000538 0.000000000347 0.000000000834 mgCcm-3soilh-1
0.0085 0.004256 0.0085 kd, d-1
0.15 0.6 0.6 kn, m3 d-1 gC-1
0.000065 0.00001625 0.000065 kl, m3 d-1 gC-1
0.000000625 0.0000025 0.0000025 kh, m3 d-1 gC-1
214000.0 19782.4 44960 Cbelowground, kg ha-1
0.00001 0.0152229 0.986 Clitter_frac, dim
0.00020 0.014 0.021 Croot_frac, dim
0.10 13710.0 29460.0 Cstem, kg ha-1
10.0 130000.0 2.0 Csapwood, kg ha-1
0.0002 0.007 0.1 Croot_coarse_frac, dim
0.002 0.000001 0.05 litter_capacity, m
0.29 0.1 0.48 theta_deep0, initial
0.28 0.05 0.48 theta_mid0, initial
0.27 0.05 0.48 theta_shallow0, initial
0.002 0.001 0.05 litter_store, initial
29.03 SLA, m2 kgC-1 leaf
1454000.0 150000.0 SRL1, m kgC-1 specific root length at root diameter of 250 um
0.004 0.064 maxRootDiam, m diameter of thickest root
0.04 minRootLifespan, years, lifespan of finest root at lowest C:N ratio
0.5 0.1 1.0 LWP_spring_minimum, -MPa
2.35 1.5 2.5 LWP_stomatal_closure, -MPa
0 is_bryophyte (1 is yes, 0 is no)
0.1 0.0 1.0 capRiseScalar, (0 to 1)
1.0 0.0 1.0 precipReduction
1.0 0.0 1.0 drainScalar, (0 to 1) proportion of drainage absorbed by water table
0.1 0.01 0.1 leafNSCscalar (proportion of leaf structural carbon)
1 usePhenology
999999999.0 leafLife Span
10 max_iteration(the_max_number_of_iterations_to_achieve_convergence_Delta<THR
7.01 0.0 2.0 microbiomeScalar, unitless, multiplier for the initial nutrient status of microbiome
0.0 snowpack_water_equivalent, m
0.0 snowpack_E_deficit_max, deg C
0.0015 melt_Rcoef, m degC-1 30-min-1
1 0 1 useHydraulics
0 0 1 useInputStress
1 0 1 useRefilling
0 0 1 forceRefilling
10 0.00001 100.0 sd_err_Ec
2.12692 0.001 15.0 sd_err_NEE
0.0 0.0 1.0 sd_err_Ec_weight

% calculate the shape for an elastic interface using surface or volume
% compressions / expansions

close all; clear

% load the parameter values

parameters_numerical;
parameters_simple;
parameters_elastic;
parameters_inverse;
rng(1); % set seed for reproducibility

% parameters for generating the artificial interface points

Nsample = 80;  % number of sample points on interface
               % NOTE: total number of points will be 2*Nsample-1
sigma_noise = 1e-4*params_phys.rneedle; % noise level for sampled points

% solve for the reference state and the deformed state

[vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
    params_num, true);

[vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
    params_num, vars_num_ref, vars_sol_ref, true);

% generate uniform data points with noise

vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
[rr_noise_ref,zz_noise_ref] = generate_noisy_shape(vars_sol_ref, ...
    vars_num_ref, Nsample, sigma_noise);

vars_sol.normals = get_normals(vars_sol, vars_num);
[rr_noise,zz_noise] = generate_noisy_shape(vars_sol, vars_num, ...
    Nsample, sigma_noise);

% fit the noisy shape with Cheby polynomials

[vars_sol_ref_fit,vars_num_ref_fit] = ...
    fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num);
vars_sol_ref_fit.p0 = vars_sol_ref.p0;

[vars_sol_fit,vars_num_fit] = ...
    fit_shape_with_chebfun(rr_noise,zz_noise,params_num);
vars_sol_fit.p0 = vars_sol.p0;

% perform CMD to find the surface stresses
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)

[vars_sol_ref_fit.sigmas, vars_sol_ref_fit.sigmap] = ...
    makeCMD(params_phys, vars_sol_ref_fit, vars_num_ref_fit);

% NOTE: in the reference state, we could also fit the YL equations to find
% the stresses. Uncomment the code below to use that approach
% [st, ~, ~, ~] = solve_inverse_young_laplace ( ...
%     vars_sol_ref_fit, params_phys, params_num, vars_num_ref_fit);
% vars_sol_ref_fit.sigmas = st*ones(vars_num.N,1);
% vars_sol_ref_fit.sigmap = st*ones(vars_num.N,1);

[vars_sol_fit.sigmas, vars_sol_fit.sigmap] = ...
    makeCMD(params_phys, vars_sol_fit, vars_num_fit);

% perform SFE to find the moduli and strains
% NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
% the the numerical results are used instead of the Cheby fit (giving a 
% best-case scenario)

[moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure,...
    vars_sol_ref_fit, vars_num_ref_fit, vars_sol_fit, vars_num_fit, ...
    params_num, true);

% post processing and plotting

errorG = abs(moduliS(1)-params_phys.Gmod)/params_phys.Gmod;
errorK = abs(moduliS(2)-params_phys.Kmod)/params_phys.Kmod;

disp(['Error in G = ', num2str(errorG*100,4), ' %']);
disp(['Error in K = ', num2str(errorK*100,4), ' %']);

plot_surface_stress(vars_num_ref.s, vars_sol_ref.sigmas, ...
    vars_sol_ref.sigmap, 2);
plot_surface_stress(vars_num_ref_fit.s, vars_sol_ref_fit.sigmas, ...
    vars_sol_ref_fit.sigmap, 2);

plot_surface_stress(vars_num.s, vars_sol.sigmas, ...
    vars_sol.sigmap, 3);
plot_surface_stress(vars_num_fit.s, vars_sol_fit.sigmas, ...
    vars_sol_fit.sigmap, 3);

plot_surface_strain(vars_num.s, vars_sol.lams, vars_sol.lamp, 4);
plot_surface_strain(vars_num_fit.s, lambda_s, lambda_r, 4);

plot_shape(rr_noise_ref, zz_noise_ref, 5);
plot_shape(vars_sol_ref_fit.r, vars_sol_ref_fit.z, 5);
plot_shape(rr_noise, zz_noise, 5);
plot_shape(vars_sol_fit.r, vars_sol_fit.z, 5);
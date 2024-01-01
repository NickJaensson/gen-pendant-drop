
close all; clear

% physical parameters
params_phys.sigma = 4;        % surface tension
params_phys.grav = 1.2;       % gravitational acceleration
params_phys.rneedle = 1.4;    % radius of the needle
params_phys.volume0 = 16;     % prescribed volume
params_phys.deltarho = 1.1;   % density difference

% numerical parameters
params_num.N = 40;          % grid points for calculation
params_num.Nplot = 80;      % grid points for plotting
params_num.eps_fw = 1e-12;  % convergence criterion forward: rms(u) < eps
params_num.maxiter = 100;   % maximum number of iteration steps

% Worthington number
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);

shape_guess = guess_shape(params_phys,1000);

vars_num = numerical_grid(params_num,[0,shape_guess.s(end)]);

vars_sol = solve_forward_young_laplace(params_phys, params_num, ...
                                       shape_guess, vars_num);

vars_num = update_numerical_grid(vars_sol,vars_num,0);

[volume,area] = calculate_volume_area(vars_sol,vars_num,1);

[s_plot,r_plot,z_plot] = interpolate_solutions(vars_sol, vars_num, ...
                                               params_num);

plot_shape(z_plot,r_plot);

[kappas,kappap] = find_curvature(vars_sol,vars_num,1);
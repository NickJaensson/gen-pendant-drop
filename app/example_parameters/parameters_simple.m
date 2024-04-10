% physical parameters for the simple droplet problem

% dimensionfull input parameters
sigma = 72;        % surface tension [mN/m]
grav = 9.807e3;    % gravitational acceleration [mm/s^2]
rneedle = 1;       % radius of the needle [mm]
volume0 = 32;      % prescribed volume in mm^3
deltarho = 1e-3;   % density difference [10^6 kg/m^3]

% dimensionless input parameters for calculation
params_phys.sigma = sigma/(deltarho*grav*rneedle^2);
params_phys.grav = 1;
params_phys.rneedle = 1;
params_phys.volume0 = volume0/rneedle^3;
params_phys.deltarho = 1;

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);
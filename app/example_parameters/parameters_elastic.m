% physical parameters for the elastic droplet problem

% dimensionfull input parameters
Kmod = 20;         % elastic dilational modulus [mN/m]
Gmod = 100;        % elastic shear modulus [mN/m]
compresstype = 1;  % 1: compress the volume other: compress the area
frac  = 0.8;       % compute elastic stresses for this compression
strainmeasure = 'pepicelli'; % which elastic constitutive model

% dimensionless input parameters for calculation
params_phys.Kmod = Kmod/(deltarho*grav*rneedle^2);
params_phys.Gmod = Gmod/(deltarho*grav*rneedle^2);
params_phys.compresstype = compresstype;  
params_phys.frac  = frac;       
params_phys.strainmeasure = strainmeasure; 
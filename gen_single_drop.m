
close all; clear

%addpath('subs/')

% physical parameters
sigma = 72;        % surface tension [mN/m]
grav = 9.807e3;    % gravitational acceleration [mm/s^2]
rneedle = 1;       % radius of the needle [mm]
volume0 = 32;      % prescribed volume in mm^3
deltarho = 1e-3;   % density difference [10^6 kg/m^3]

% numerical parameters
N = 40;          % resolution of the discretization for calculation
Nplot = 80;      % resolution of the discretization for plotting
Ncheb = 10;      % number of Chebyshev to describe the shape
alpha = 1;     % relaxation parameter in the Newton-Raphson scheme

tic

% NOTE: the calculation is done in dimensionless form, using the 
% dimensionless surface tension sigma' and volume V'

% calculate the dimensionless quantities
sigmaprime = sigma/(deltarho*grav*rneedle^2);
volume0prime = volume0/rneedle^3;

% find the initial guess of the droplet shape
if  deltarho*grav*volume0/(2*pi*sigma*rneedle) > 0.14

    % predict the droplet shape using the emperical approach from Nagel
    
    % predict the maximum length of the interface (empirical Nagel)
    smax = sqrt(sigmaprime)*2.0/0.8701;

    % get the differentation/integration matrices and the grid
    [D,~,w,s] = dif1D('cheb',0,smax,N,5);

    % predict the shape of the interface (empirical Nagel)
    z = -4/3*smax/pi*(cos(pi*3/4*s/smax));
    z = z - max(z);
    r = 4/3*smax/pi*(sin(pi*3/4*s/smax));
    psi = pi*3/4*s/smax;

    C = 1; % initial stretch parameter
    p0 = sqrt(sigmaprime)*1.5; % predict the pressure (empirical Nagel)

else
    
    % predict the droplet shape using a quarter of a period of a cosine
    % with similar volume as imposed

    % find the initial guess of the droplet shape
    Ntemp = 1000;
    Rtemp = 1.0; 
    r_guess = linspace(0,Rtemp,Ntemp);
    z_guess = -sqrt(2*volume0prime/(pi*Rtemp))*cos(pi*r_guess/(2*Rtemp));

    % determine the curve length
    ds_guess = zeros(size(r_guess));
    ds_guess(2:end) = sqrt((r_guess(2:end)-r_guess(1:Ntemp-1)).^2 + ...
                           (z_guess(2:end)-z_guess(1:Ntemp-1)).^2);
    s_guess = cumsum(ds_guess);
    smax = s_guess(end);

    % get the differentation/integration matrices and the grid
    [D,~,w,s] = dif1D('cheb',0,smax,N,5);

    % interpolate the shape in the Chebyshev points
    r = interp1(s_guess,r_guess,s);
    z = interp1(s_guess,z_guess,s);

    psi = atan2(D*z,D*r);      % intial psi value 
    C = 1;                     % initial stretch parameter
    p0 = 2*Rtemp*sigmaprime;   % predict the pressure
    
end
    
% initialize some variables 
Z = zeros(N);            % matrix filled with zeros
IDL = [1, zeros(1,N-1)]; % line with single one and rest zeros
ZL = zeros(1,N);         % line completely filled with zeros
u = ones(3*N+2,1); b = ones(3*N+2,1); % solution vector and right hand side
iter = 0; crash = 0; 

while rms(u) > 1e-10

  iter = iter + 1;
  
  if iter > 1200 
    warning('iter > 12000!');
    crash = 1; break;
  end

  % determine r from psi
  A11 = C*D; A13 = diag(sin(psi)); A14 = D*r; b1 = -(C*D*r-cos(psi));

  % determine z from psi 
  A22 = C*D; A23 = diag(-cos(psi)); A24 = D*z; b2 = -(C*D*z-sin(psi));

  % determine psi from Laplace law
  A31 = -sigmaprime*diag(sin(psi)./r.^2);
  A32 = diag(ones(N,1));
  A33 = C*sigmaprime*D + sigmaprime*diag(cos(psi)./r);
  A34 = sigmaprime*(D*psi);
  A35 = -ones(N,1);
  b3 = p0-z-sigmaprime*(C*D*psi+sin(psi)./r);

  % impose the needle radius as a BC (imposes the domain length)
  % NOTE: the lengths are scaled with the radius, thus its value is one
  A41 = fliplr(IDL); b4 = (1-r(end));
  
  % determine pressure - use volume
  A51 = pi*2*w.*r'.*sin(psi');
  A53 = pi*w.*r'.^2.*cos(psi');
  A54 = -volume0prime;
  b5 = -(pi*w*(r.^2.*sin(psi))-C*volume0prime);

  % boundary condition r(0) = 0
  A11(1,:) = IDL; 
  A13(1,:) = ZL; 
  A14(1) = 0;
  b1(1) = -r(1);
  
  % boundary condition z(s0) = 0
  A22(1,:) = fliplr(IDL); 
  A23(1,:) = ZL; 
  A24(1) = 0;
  b2(1) = -z(end);
  
  % boundary condition phi(0) = 0
  A31(1,:) = ZL; 
  A32(1,:) = ZL; 
  A33(1,:) = IDL; 
  A34(1,:) = 0; 
  A35(1,:) = 0;
  b3(1) = -psi(1);

  % assemble matrices
  Z1 = zeros(N,1);
     
  A = [[A11,   Z, A13, A14,  Z1];
       [  Z, A22, A23, A24,  Z1];
       [A31, A32, A33, A34, A35];
       [A41,  ZL,  ZL,   0,   0];
       [A51, Z1', A53, A54,   0]];
     
  b = [b1;b2;b3;b4;b5];

  % solve the system of equations
  u = A\b;

  % update variables
  r   = r   + alpha*u(1:N);
  z   = z   + alpha*u(N+1:2*N);
  psi = psi + alpha*u(2*N+1:3*N); 
  C   = C   + alpha*u(3*N+1);
  p0  = p0  + alpha*u(3*N+2);

  if rms(b) > 1e3
    crash = 1; break;
  end

  fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

end

toc

% compute volume and area (scaled back to dimensionfull)
disp(['volume = ', num2str(rneedle^3*pi*w*(r.^2.*sin(psi))/C,15),' mm^3']);
disp(['area = ', num2str(rneedle^2*pi*2*w*(r)/C,15),' mm^2']);
disp(['pressure = ', num2str(deltarho*grav*rneedle*p0,15),' Pa']);

% % plot the shape of the drop on the numerical grid
% figure; hold on
% scatter(rneedle*r',rneedle*z','b');
% plot(rneedle*r',rneedle*z','b');
% set(gca,'DataAspectRatio',[1 1 1])

% interpolate the numerical solutions on a finer grid. 
% NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
% though all the points (see book of Trefethen on Spectral Methods in 
% Matlab, page  63). For plotting purposes we use a simpler interpolation 
ss = linspace(s(1),s(end),Nplot)';
rr = interp1(s,r,ss,'pchip');
zz = interp1(s,z,ss,'pchip');

% plot the shape of the drop on the plotting grid
figure; hold on
scatter(rneedle*rr',rneedle*zz','b');
plot(rneedle*rr',rneedle*zz','b');
set(gca,'DataAspectRatio',[1 1 1])
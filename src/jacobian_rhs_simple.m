function [A,b] = jacobian_rhs_simple(params,itervars)
    
    D = params.D;
    w = params.w;
    r = itervars.r;
    z = itervars.z;
    psi = itervars.psi;
    C = itervars.C;
    p0 = itervars.p0;
    
    % initialize some variables 
    Z = zeros(params.N);            % matrix filled with zeros
    IDL = [1, zeros(1,params.N-1)]; % line with single one and rest zeros
    ZL = zeros(1,params.N);         % line completely filled with zeros 
    b = ones(3*params.N+2,1); % solution vector and right hand side
    
    % determine r from psi
    A11 = C*D; A13 = diag(sin(psi)); A14 = D*r; b1 = -(C*D*r-cos(psi));
    
    % determine z from psi 
    A22 = C*D; A23 = diag(-cos(psi)); A24 = D*z; b2 = -(C*D*z-sin(psi));
    
    % determine psi from Laplace law
    A31 = -params.sigma*diag(sin(psi)./r.^2);
    A32 = diag(ones(params.N,1));
    A33 = C*params.sigma*D + params.sigma*diag(cos(psi)./r);
    A34 = params.sigma*(D*psi);
    A35 = -ones(params.N,1);
    b3 = p0-z-params.sigma*(C*D*psi+sin(psi)./r);
    
    % impose the needle radius as a BC (imposes the domain length)
    % NOTE: the lengths are scaled with the radius, thus its value is one
    A41 = fliplr(IDL); b4 = (params.rneedle-r(end));
    
    % determine pressure - use volume
    A51 = pi*2*w.*r'.*sin(psi');
    A53 = pi*w.*r'.^2.*cos(psi');
    A54 = -params.volume0;
    b5 = -(pi*w*(r.^2.*sin(psi))-C*params.volume0);
    
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
    Z1 = zeros(params.N,1);
     
    A = [[A11,   Z, A13, A14,  Z1];
       [  Z, A22, A23, A24,  Z1];
       [A31, A32, A33, A34, A35];
       [A41,  ZL,  ZL,   0,   0];
       [A51, Z1', A53, A54,   0]];
     
    b = [b1;b2;b3;b4;b5];


end
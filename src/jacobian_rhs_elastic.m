function [A, b] = jacobian_rhs_elastic(params_phys, vars_sol, vars_num)
    % JACOBIAN_RHS_ELASTIC Computes the Jacobian matrix and RHS vector for 
    % the elastic shape problem.
    %
    % INPUTS:
    %   params_phys - Structure with physical parameters
    %   vars_sol    - Structure with solution variables
    %   vars_num    - Structure with numerical variables
    %
    % OUTPUTS:
    %   A - Jacobian matrix for the elastic system
    %   b - Right-hand side vector for the elastic system

    % local copy of diffmat and intmat
    % NOTE: D and w are wrt to the numerical reference domain (described by
    % s0, but we use D and w here for convenience
    D = vars_num.D0;
    w = vars_num.w0;

    % local copy of other variables
    r = vars_sol.r;
    z = vars_sol.z;
    psi = vars_sol.psi;
    sigmas = vars_sol.sigmas;
    sigmap = vars_sol.sigmap;
    lams = vars_sol.lams;
    lamp = vars_sol.lamp;
    p0 = vars_sol.p0;
    N = vars_num.N;
    C = vars_num.C;

    Kmod = params_phys.Kmod;
    Gmod = params_phys.Gmod;

    % initialize some variables 
    Z = zeros(N);            % matrix filled with zeros
    IDL = [1, zeros(1,N-1)]; % line with single one and rest zeros
    ZL = zeros(1,N);         % line completely filled with zeros
    Z1 = zeros(N,1);         % column filled with zeros
    I = eye(N);              % unit matrix

    % Eq. 1-4 Knoche, p85, eq.5.7, Eq. 5-7, eq.5.8
    % determine r from psi (incl lams)
    % A11 = C*D
    % A13 = lams*sin(psi)
    % A16 = -(C*D*r)/lams
    % b1 = lams*cos(psi) - (C*D*r)
    A11 = C*diag(1./lams)*D;
    A13 = diag(sin(psi));
    A16 = -C*diag((D*r)./(lams.^2));
    b1 = cos(psi)-C*(D*r)./lams;

    % determine z from psi (incl lams)
    % A22 = C*D
    % A23 = -lams*cos(psi)
    % A26 = -(C*D*z)/lams
    % b2 = lams*sin(psi) - (C*D*z)
    A22 = C*diag(1./lams)*D;
    A23 = diag(-cos(psi));
    A26 = -C*diag((D*z)./(lams.^2));
    b2 = sin(psi)-C*(D*z)./lams;

    % determine psi from laplace law
    % A31 = -lams*(sigmap*sin(psi))/r^2
    % A32 = lams*g*rho
    % A33 = (lams*sigmap*cos(psi))/r + (C*sigmas)*D
    % A34 = (C*D*psi)
    % A35 = lams*sin(psi)/r
    % A36 = -(C*D*psi*sigmas)/lams
    % A38 = -lams
    % b3 = lams*P - lams*(sigmap*sin(psi))/r - lams*g*rho*z - (C*D*psi*sigmas)
    A31 = -diag(sigmap.*sin(psi)./(r.^2));
    A32 = eye(N)*params_phys.deltarho*params_phys.grav;
    A33 = diag(sigmap.*cos(psi)./r)+C*diag(sigmas./lams)*D;
    A34 = C*diag((D*psi)./lams);
    A35 = diag(sin(psi)./r);
    A36 = -C*diag((D*psi).*sigmas./(lams.^2));
    A38 = -ones(N,1);
    b3 = p0 - sigmap.*sin(psi)./r ...
                  - z*params_phys.deltarho*params_phys.grav -C*sigmas.*(D*psi)./lams;


    if params_phys.compresstype == 1
        % determine pressure - use volume
        % A81 = (2*int*lams*r*pi*sin(psi))
        % A83 = (int*lams*r^2*pi*cos(psi))
        % A86 = (int*r^2*pi*sin(psi))
        % b8 = V*C - (int*lams*r^2*pi*sin(psi))        
        wdef = w.*lams'/C; 
        A81 = 2*pi*wdef.*(r.*sin(psi))';
        A83 =   pi*wdef.*(r.^2.*cos(psi))';
        A86 =   pi*(w/C).*((r.^2).*sin(psi))';
        b8 =   -pi*wdef*((r.^2).*sin(psi))+params_phys.volume;
    else
        % determine pressure - use area
        % A81 = (2*int*lams*pi)/C
        % A86 = (2*int*r*pi)/C
        % b8 = A - (2*int*lams*r*pi)/C
        wdef = w.*lams'/C; 
        A81 = 2*pi*wdef;
        A83 =   zeros(1,N);
        A86 =   2*pi*(w/C).*r';
        b8 =   -2*pi*wdef*r+params_phys.area;
    end

    % boundary condition r(0) = 0
    A11(1,:) = IDL;
    A13(1,:) = ZL;
    A16(1,:) = ZL;
    b1(1) = -r(1);

    % boundary condition z(s0) = 0
    A22(1,:) = fliplr(IDL);
    A23(1,:) = ZL;
    A26(1,:) = ZL;
    b2(1) = -z(end);

    % boundary condition phi(0) = 0
    A31(1,:) = ZL;
    A32(1,:) = ZL;
    A33(1,:) = IDL;
    A34(1,:) = ZL;
    A35(1,:) = ZL;
    A36(1,:) = ZL;
    A38(1,:) = 0;
    b3(1) = -psi(1);

    % determine sigmas from projection of force balance
    % A41 = (C*D*sigmas)
    % A43 = -lams*sin(psi)*(sigmas - sigmap)
    % A44 = lams*cos(psi) + (C*r)*D
    % A45 = -lams*cos(psi)
    % A46 = -(C*r*D*sigmas)/lams
    % b4 = - lams*cos(psi)*(sigmas - sigmap) - (C*r*D*sigmas)    
    A41 = C*diag((D*sigmas)./lams);
    A43 = diag(sin(psi).*(sigmap-sigmas));
    A44 = diag(cos(psi))+C*diag(r./lams)*D;
    A45 = -diag(cos(psi));
    A46 = -C*diag(r.*(D*sigmas)./(lams.^2));
    b4 = -C*r.*(D*sigmas)./lams+cos(psi).*(sigmap-sigmas);

    % define some convenient variables for constitutive equation
    lamsm1 = lams.^(-1); lamsm2 = lams.^(-2); lamsm3 = lams.^(-3);
    lampm1 = lamp.^(-1); lampm2 = lamp.^(-2); lampm3 = lamp.^(-3);
    J = lams.*lamp;
    K = params_phys.Kmod;
    G = params_phys.Gmod;
        
    switch params_phys.strainmeasure

    case 'linear_hookean'

        A54 = eye(N);
        A56 = (-G-K)*eye(N);
        A57 =  (G-K)*eye(N);
        b5 = params_phys.sigma - sigmas + (lams - 1)*(G + K) - (G - K)*(lamp - 1);

        A65 = eye(N);
        A66 = (G-K)*eye(N);
        A67 = (-G-K)*eye(N);
        b6 = params_phys.sigma - sigmap + (lamp - 1)*(G + K) - (G - K)*(lams - 1);
        
    case 'hencky'

        A54 = eye(N);
        A56 = diag(-G*lamsm1 -K*lamsm1);
        A57 = diag( G*lampm1 -K*lampm1);
        b5 = params_phys.sigma - sigmas + K*log(J) + G*log(lams.*lampm1);

        A65 = eye(N);
        A66 = diag( G*lamsm1 -K*lamsm1);
        A67 = diag(-G*lampm1 -K*lampm1);
        b6 = params_phys.sigma - sigmap + K*log(J) + G*log(lamp.*lamsm1);

    case 'pepicelli'

        % determine sigma^p
        A54 = eye(N);
        A56 = diag(K*log(J).*lamsm2.*lampm1 - K*lamsm2.*lampm1 - G*lamsm3);
        A57 = diag(G*lampm3 - K*lamsm1.*lampm2 + K*log(J).*lamsm1.*lampm2);
        b5 = params_phys.sigma - sigmas - G*(lamsm2-lampm2)/2 + K*log(J)./J;
        
        % determine lambda^s
        A65 = eye(N);
        A66 = diag(K*log(J).*lamsm2.*lampm1 - K*lamsm2.*lampm1 + G*lamsm3);
        A67 = diag(-G*lampm3 - K*lamsm1.*lampm2 + K*log(J).*lamsm1.*lampm2);
        b6 = params_phys.sigma - sigmap + G*(lamsm2-lampm2)/2 + K*log(J)./J;

    end

    % A71 = 1
    % A77 = -rstar
    % b7 = -r + lamp*rstar
    % determine lambda^r
    A71 = eye(N);
    A77 = -diag(vars_sol.r_star);
    b7 = -r+lamp.*vars_sol.r_star;

    % boundary condition dsigmas/ds(0) = 0
    % NOTE: this BC is included in the Newton-Raphson iteration
    A41(1,:) = ZL;
    A43(1,:) = ZL;
    A44(1,:) = (1/lams(1))*vars_num.C*D(1,:);
    A45(1,:) = ZL;
    A46(1,:) = ZL;
    A46(1,1) = -(1/lams(1)^2)*vars_num.C*(D(1,:)*sigmas);
    b4(1) = -(1/lams(1))*vars_num.C*(D(1,:)*sigmas);

    % boundary condition lamp(s0) = 1
    A71(end,:) = ZL;
    A77(end,:) = fliplr(IDL);
    b7(end) =  1.0 - lamp(end);

    % boundary condition sigmas(0)=sigmap(0)
    A54(1,:) = IDL;
    A56(1,:) = 0;
    A57(1,:) = 0;
    A55 = Z;
    A55(1,:) = -IDL;
    b5(1) = 0;

    % boundary condition lams(0)=lamp(0)
    A71(1,:) = 0;
    A77(1,:) = IDL;
    A76 = Z;
    A76(1,:) = -IDL;
    b7(1) = 0;

    % combine matrices
    A = [[A11,   Z, A13,   Z,    Z, A16,   Z,  Z1]; ...
         [  Z, A22, A23,   Z,    Z, A26,   Z,  Z1]; ...
         [A31, A32, A33, A34,  A35, A36,   Z, A38]; ...
         [A41,   Z, A43, A44,  A45, A46,   Z,  Z1]; ...
         [  Z,   Z,   Z, A54,  A55, A56, A57,  Z1]; ...
         [  Z,   Z,   Z,   Z,  A65, A66, A67,  Z1]; ...
         [A71,   Z,   Z,   Z,    Z, A76, A77,  Z1]; ...
         [A81, Z1', A83,  Z1', Z1', A86, Z1',   0]];

    b = [b1;b2;b3;b4;b5;b6;b7;b8];

end

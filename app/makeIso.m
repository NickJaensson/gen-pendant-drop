function [ tension, pcap, rrlaplace, zzlaplace ] = makeIso(zz_in,rr_in,psi_in,diffmat_in)
    % makeIso(toplot) = [tension, pstat] fits the shape functions to surface 
    % tension and pressure. RMS fitting of R using the Schur complement.
    % toplot specifies if the result is to be plotted.

    % assign some local variables
    psi = psi_in;
    r = rr_in;
    z = zz_in;
    d = diffmat_in;

    g_useP = 0; p_measured = 0;

    N = length(r);
    C = 1;
    Gam = 10;
    if g_useP
        P = p_measured;
    else
        P = Gam/2;
    end
    alpha = 0.25;
    iter = 1;u=1;
    rmsu = 1e3;

    while rms(u) > 1e-9
        
        iter = iter+1;
        
        if iter > 100
            error('Iteration did not converge!')
        end  

        % create matrix
        [tA, tb] =  matrix_iso(0,P,Gam,d,0,C,r,z,psi);

        % pressure solved or prescribed
        if g_useP
            A2 = tA(:,end-1);
        else
            A2 = tA(:,end-1:end);
        end

         % Schur Element
        A1 = tA(:,1:end-2);

        A3 = [diag(rr_in-r), zeros(1*N,2*N)];

        IA1 = inv(A1);
        SchurA = A3*(IA1*A2);
        b = A3*(IA1*tb)-(rr_in-r).^2;
        
        if ( iter==2 ); warning('off', 'MATLAB:rankDeficientMatrix'); end
        u2 = SchurA\b;
        if ( iter==2 ); warning('on', 'MATLAB:rankDeficientMatrix'); end

        u1 = IA1*(tb-A2*u2);
        u = [u1;u2];

        % update variables
        r = r+alpha*u(1:N);
        z = z+alpha*u(N+1:2*N);
        psi = psi+alpha*u(2*N+1:3*N);
        Gam = Gam+alpha*u(3*N+1);
        if ~g_useP
            P = P+alpha*u(end);
        end

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

        rmsb = rms(tb);
    end

    rrlaplace = r;
    zzlaplace = z;
    tension = Gam;
    pcap = P;
    
end

function [ A, b] = matrix_iso(~, P, Gam, d,~,C,r,z,psi)

    % matrix and rhs for isotropic interface,
    % its unknowns are: Gam, P.

    % Final check uses this routine in all functions for drop create and detect

    N = length(r);

    % full problem in r,z and psi
    % zero matrix, zero line and identity line array
    Z = zeros(N);
    IDL = [1, zeros(1,N-1)];
    ZL = zeros(1,N);

    % determine r from psi
    A11 = C*d; % N x N
    A13 = diag(sin(psi)); % N x N
    A14 = [zeros(N,1), zeros(N,2)]; % N x 3
    b1 = cos(psi)-C*d*r;  % N x 1

    % boundary condition r(1) = 0
    A11(1,:) = IDL; 
    A13(1,:) = ZL;
    A14(1,:) = zeros(1,3);
    b1(1) = -r(1);

    % determine z from psi
    A22 = C*d;  % N x N
    A23 = diag(-cos(psi)); % N x N
    A24 = [zeros(N,1), zeros(N,2)]; % N x 3
    b2 = sin(psi)-C*d*z; % N x 1

    % boundary condition z(end) =0
    A22(end,:) = fliplr(IDL);
    A23(end,:) = ZL;
    A24(end,:) = zeros(1,3);
    b2(end) = -z(end);

    % determine psi from Laplace law
    A31 = diag(-sin(psi)./r.^2); % N x N
    A32 = eye(N); % N x N
    A33 = C*Gam*d+diag(Gam*cos(psi)./r); % N x N
    A34 = [Gam*d*psi,d*psi+sin(psi)./r, -ones(N,1) , zeros(N,0)]; % N x 3
    b3 = -z+P-C*Gam*(d*psi)-Gam*sin(psi)./r; % N x 1

    % boundary condition phi(0) = 0
    A31(1,:) = ZL;
    A32(1,:) = ZL;
    A33(1,:) = IDL;
    A34(1,:) = zeros(1,3);
    b3(1) = -psi(1);

    % Build small matrix, recheck.

    A = [[A11, Z, A13, A14(:,2:end)];[Z, A22, A23, A24(:,2:end)];[A31, A32, A33, A34(:,2:end)]];
    b = [b1;b2;b3];
    
end
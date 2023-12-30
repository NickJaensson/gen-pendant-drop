function [vars_sol,params] = solve_forward_young_laplace(params)

    % calculate the Worthinton number
    params.Wo = params.deltarho*params.grav*params.volume0/...
                                        (2*pi*params.sigma*params.rneedle);

    % find an initial guess of the shape
    [r_guess, z_guess, s_guess] = guess_shape(params,1000);
    
    smax = s_guess(end); % the total length of the 1D domain
    
    % get the differentation/integration matrices and the grid
    [params.D,~,params.w,params.s] = numerical_grid(params.N,[0,smax]);
    
    % interpolate the shape in the Chebyshev points
    r = interp1(s_guess,r_guess,params.s);
    z = interp1(s_guess,z_guess,params.s);
    
    psi = atan2(params.D*z,params.D*r);   % intial psi value 
    C = 1;                                % initial stretch parameter
    p0 = 2*params.sigma/params.rneedle;   % initial pressure
    u = ones(3*params.N+2,1);             % initial solution vector
    
    % store some variables for the iteration
    iter = 0;
    vars_sol.r = r; vars_sol.z = z; vars_sol.psi = psi;
    vars_sol.C = C; vars_sol.p0 = p0; 
    
    % start the Newton-Raphson iteration
    while rms(u) > params.eps
    
        iter = iter + 1;
        
        if iter > params.maxiter
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_simple(params,vars_sol);
        
        % solve the system of equations
        u = A\b;
        
        % update variables
        vars_sol.r   = vars_sol.r   + params.alpha*u(1:params.N);
        vars_sol.z   = vars_sol.z   + params.alpha*u(params.N+1:2*params.N);
        vars_sol.psi = vars_sol.psi + params.alpha*u(2*params.N+1:3*params.N); 
        vars_sol.C   = vars_sol.C   + params.alpha*u(3*params.N+1);
        vars_sol.p0  = vars_sol.p0  + params.alpha*u(3*params.N+2);    
        
        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));
    
    end

    % the integration and differentation matrices in the solution state
    params.ws = params.w/vars_sol.C; 
    params.Ds = vars_sol.C*params.D; 
    
end
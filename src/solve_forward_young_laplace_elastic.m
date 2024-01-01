function [vars_sol, vars_num] = solve_forward_young_laplace_elastic(vars_sol, params_phys, params_num, vars_num)
% SOLVE_FORWARD_YOUNG_LAPLACE_ELASTIC Solves the elastic Young-Laplace 
% equation for the shape of a drop, for different compressions
%
% This function solves the elastic Young-Laplace equation for a constant 
% surface tension and elastic moduli, using a Newton-Raphson iteration. 
% Elastic stresses are generated by imposing volume compressions
%
% Usage:
%   [vars_sol, vars_num] = solve_forward_young_laplace_elastic(vars_sol, params_phys, params_num, vars_num)
%
% Inputs:
%   vars_sol - Initial structure containing variables of the solution state.
%   params_phys - A structure containing physical parameters:
%                 volume0 (initial volume), frac (volume fraction for compression),
%                 compresstype (type of compression, currently only '1' for volume compression is implemented).
%   params_num - A structure containing numerical parameters:
%                N (number of points), eps_fw (tolerance for convergence),
%                maxiter (maximum iterations), alpha (relaxation parameter).
%   vars_num - A structure containing numerical variables used in the solution of the reference state:
%              Differentiation and integration matrices, Chebyshev points, and modified
%              matrices post-solution (ws, Ds).
%
% Outputs:
%   vars_sol - Updated structure containing the solution variables after convergence:
%              r, z (coordinates of the interface), psi (interface angle),
%              sigmas, sigmap (stress components), lams, lamp (deformation components),
%              p0 (pressure).
%   vars_num - Updated structure containing numerical variables used in the solution:
%              Differentiation and integration matrices, Chebyshev points, and modified
%              matrices post-solution (sdef, wdef, Ddef).

    % determine the current target volume/area
    if params_phys.compresstype == 1
        params_phys.volume = params_phys.volume0*params_phys.frac;
    else
        error('area compression not implemented')    
    end

    % store some variables for the iteration
    iter = 0; u = ones(3*params_num.N+2,1);

    % start the Newton-Raphson iteration
    while rms(u) > params_num.eps_fw
    
        iter = iter + 1;
        
        if iter > params_num.maxiter
            error('Iteration did not converge!')
        end    
    
        % build the Jacobian and RHS
        [A,b] = jacobian_rhs_elastic(params_phys,vars_sol,vars_num);
        
        % solve the system of equations
        u = A\b;
    
        % update variables
        vars_sol.r   = vars_sol.r + params_num.alpha*u(1:params_num.N);
        vars_sol.z   = vars_sol.z + params_num.alpha*u(params_num.N+1:2*params_num.N);
        vars_sol.psi = vars_sol.psi + params_num.alpha*u(2*params_num.N+1:3*params_num.N);    
        vars_sol.sigmas = vars_sol.sigmas + params_num.alpha*u(3*params_num.N+1:4*params_num.N);
        vars_sol.sigmap = vars_sol.sigmap + params_num.alpha*u(4*params_num.N+1:5*params_num.N);
        vars_sol.lams = vars_sol.lams + params_num.alpha*u(5*params_num.N+1:6*params_num.N);
        vars_sol.lamp = vars_sol.lamp + params_num.alpha*u(6*params_num.N+1:7*params_num.N);
        vars_sol.p0  = vars_sol.p0 + params_num.alpha*u(end);

        fprintf('iter %d: rms(u) = %d\n',iter,rms(u));

    end

    % the integration and differentation matrices in the deformed state
    % NOTE: this construction of Ddef is simlar to first applying D*f/C,
    % and then dividing the components by the components of (1/lams)
    vars_num.wdef = vars_num.w.*vars_sol.lams'/vars_num.C; 
    vars_num.Ddef = vars_num.C*vars_num.D.*repelem((1./vars_sol.lams)',params_num.N,1); 

    % construct the integration matrix from the integration vector
    vars_num.wmat = repmat(vars_num.w,params_num.N,1);
    vars_num.wmat = tril(vars_num.wmat);

    % compute the value of s in the deformed state
    vars_num.sdef = vars_num.wmat*vars_sol.lams/vars_num.C;
   
end
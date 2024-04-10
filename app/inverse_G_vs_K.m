% calculate the Laplace shape for a given surface tension and given
% pressure/volume/area

close all; clear

% load the parameter values

% numerical parameters
params_num.N = 40;          % grid points for calculation
params_num.eps_fw_simple = 1e-12;  % convergence criterion forward
params_num.maxiter_simple = 100;   % maximum number of iteration steps
params_num.eps_fw_elastic = 1e-10;  % convergence criterion forward
params_num.maxiter_elastic = 100;   % maximum number of iteration steps

% numerical parameters for inverse problem
% NOTE: it seems that this must be increased for droplets that are more
% spherical (e.g. 1e-3 for Vo=5 and 1e-2 for Vo=1)
params_num.eps_cheb = 1e-3;      % error for describing the shape

% parameters for SFE inverse problem
params_num.eps_inv_sfe = 1e-6;   % convergence critertion inverse problem
params_num.K_guess = 20;          % guess for dilational modulus
params_num.G_guess = 20;          % guess for shear moduls
params_num.alpha_sfe = 0.1;      % relaxation parameter in inverse problem
params_num.maxiter_inv_sfe = 5000;% maximum number of iteration steps invers
 
rng(1); % set seed for reproducibility

Vo = 5;  % called Ar in the manuscript
Bo = 1;  % called Wo in the manuscript

% physical parameters for the simple droplet problem
params_phys.sigma = 60;       % surface tension
params_phys.grav = 9.807e3;   % gravitational acceleration
params_phys.deltarho = 1e-3;  % density difference

params_phys.rneedle = sqrt(Bo * params_phys.sigma / ...
    (params_phys.deltarho * params_phys.grav * Vo));
params_phys.volume0 = Vo*params_phys.rneedle^3;

% Worthington number (needed for initial shape guess)
params_phys.Wo = params_phys.deltarho*params_phys.grav*...
    params_phys.volume0/(2*pi*params_phys.sigma*params_phys.rneedle);

% physical parameters for the elastic problem
params_phys.compresstype = 2;  % 1: compress the volume other: compress the area
params_phys.frac  = 0.9;       % compute elastic stresses for this compression
params_phys.strainmeasure = 'pepicelli'; % which elastic constitutive model

K_all = 1:10:100;    % dilational modulus
G_all = 1:10:100;    % shear modulus

% parameters for generating the artificial interface points

Nsample = 80;  % number of sample points on interface
               % NOTE: total number of points will be 2*Nsample-1
sigma_noise = 1e-4*params_phys.rneedle; % noise level for sampled points

error_mean_G = zeros(length(K_all),length(G_all));
error_mean_K = zeros(length(K_all),length(G_all));
nr_of_samples = zeros(length(K_all),length(G_all));

for iii = 1:length(K_all)
    
    for jjj = 1:length(G_all)

        Nrand = 100; % nunber of random samples

        params_phys.Kmod = K_all(iii);
        params_phys.Gmod = G_all(jjj);
        
        % solve for the reference state and the deformed state
        
        [vars_num_ref, vars_sol_ref, params_phys] = gen_single_drop(params_phys, ...
            params_num, false);
        
        % initialize current inferred values

        G{iii,jjj} = [NaN];
        K{iii,jjj} = [NaN];

        try
            [vars_num, vars_sol] = gen_single_drop_elastic(params_phys, ...
                params_num, vars_num_ref, vars_sol_ref, false);
            vars_sol_ref.normals = get_normals(vars_sol_ref, vars_num_ref);
            vars_sol.normals = get_normals(vars_sol, vars_num);
            drop_succesful = 1;
        catch
            drop_succesful = 0;
            warning('forward problem failed');
            fprintf('Gmod = %d  Kmod = %d\n',[params_phys.Gmod,params_phys.Kmod]);
        end           


        if drop_succesful

            kk = 1;
    
            for k = 1:Nrand
    
                try
    
                    % generate uniform data points with noise
                    
                    [rr_noise_ref,zz_noise_ref] = generate_noisy_shape(vars_sol_ref, ...
                        vars_num_ref, Nsample, sigma_noise);
                    
                    [rr_noise,zz_noise] = generate_noisy_shape(vars_sol, vars_num, ...
                        Nsample, sigma_noise);
    
                    % fit the noisy shape with Cheby polynomials
                    
                    [vars_sol_ref_fit{iii,jjj}{kk},vars_num_ref_fit] = ...
                        fit_shape_with_chebfun(rr_noise_ref,zz_noise_ref,params_num);
                    vars_sol_ref_fit{iii,jjj}{kk}.p0 = vars_sol_ref.p0;
                    
                    [vars_sol_fit{iii,jjj}{kk},vars_num_fit] = ...
                        fit_shape_with_chebfun(rr_noise,zz_noise,params_num);
                    vars_sol_fit{iii,jjj}{kk}.p0 = vars_sol.p0;
                    
                    % perform CMD to find the surface stresses
                    % NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
                    % the the numerical results are used instead of the Cheby fit (giving a 
                    % best-case scenario)
                    
                    [vars_sol_ref_fit{iii,jjj}{kk}.sigmas, vars_sol_ref_fit{iii,jjj}{kk}.sigmap] = ...
                        makeCMD(params_phys, vars_sol_ref_fit{iii,jjj}{kk}, vars_num_ref_fit);
                    
                    [vars_sol_fit{iii,jjj}{kk}.sigmas, vars_sol_fit{iii,jjj}{kk}.sigmap] = ...
                        makeCMD(params_phys, vars_sol_fit{iii,jjj}{kk}, vars_num_fit);
    
                    % perform SFE to find the moduli and strains
                    % NOTE: by replacing vars_sol_fit -> vars_sol and vars_num_fit -> vars_num
                    % the the numerical results are used instead of the Cheby fit (giving a 
                    % best-case scenario)
                    
                    [moduliS, lambda_s, lambda_r]  = makeSFE(params_phys.strainmeasure,...
                        vars_sol_ref_fit{iii,jjj}{kk}, vars_num_ref_fit, vars_sol_fit{iii,jjj}{kk}, vars_num_fit, ...
                        params_num, false);
                    
                    % post processing and plotting
                    
                    G{iii,jjj}(kk) = moduliS(1);
                    K{iii,jjj}(kk) = moduliS(2);
                    kk = kk + 1;

                    if kk == 11
                        break
                    end
    
                catch
    
                end
    
            end
        end

        G_error = abs((G{iii,jjj}-params_phys.Gmod)/params_phys.Gmod);
        K_error = abs((K{iii,jjj}-params_phys.Kmod)/params_phys.Kmod);

        error_mean_G(iii,jjj) = mean(G_error);
        error_mean_K(iii,jjj) = mean(K_error);
        nr_of_samples(iii,jjj) = kk;

    end
end

% plot the number of samples for each G and K
imagesc(K_all, G_all, nr_of_samples);
colorbar
set(gca,'YDir','normal');
set(gca, 'XTick', K_all, 'XTickLabel', K_all);
set(gca, 'YTick', G_all, 'YTickLabel', G_all);

% plot the errors with the NaNs removed
figure
imAlpha = ones(size(error_mean_G));
imAlpha(isnan(error_mean_G))=0;
imagesc(K_all, G_all, error_mean_G,'AlphaData',imAlpha);
colorbar
set(gca,'color',1.0*[1,1,1]);
set(gca,'ColorScale','log');
set(gca,'YDir','normal');
set(gca, 'XTick', K_all, 'XTickLabel', K_all);
set(gca, 'YTick', G_all, 'YTickLabel', G_all);

% plot the errors with the NaNs removed
figure
imAlpha = ones(size(error_mean_K));
imAlpha(isnan(error_mean_K))=0;
imagesc(K_all, G_all, error_mean_K,'AlphaData',imAlpha);
colorbar
set(gca,'color',1.0*[1,1,1]);
set(gca,'ColorScale','log');
set(gca,'YDir','normal');
set(gca, 'XTick', K_all, 'XTickLabel', K_all);
set(gca, 'YTick', G_all, 'YTickLabel', G_all);
close all; clear

gen_single_drop_elastic;

close all

% disp ( abs(volume-12.8000000000145) );
% disp ( abs(area-22.5156483902244) );
% disp ( abs(p0-2.02056164104927) );
% disp ( abs(max(sigmas)-3.33958761227839) );
% disp ( abs(max(sigmap)-3.86864491619739) );

% compare to old values (gen-pendant-drop before refactoring:
eps_test = 1e-10; 
assert ( abs(volume-12.8000000000001) < eps_test );
assert ( abs(area-22.5156483902096) < eps_test );
assert ( abs(vars_sol.p0-2.02056164104124) < eps_test );
assert ( abs(max(vars_sol.sigmas)-3.33958761227925) < eps_test );
assert ( abs(max(vars_sol.sigmap)-3.86864491619756) < eps_test );

disp('All tests passed!')
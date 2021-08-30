close all
clear

%load eigenpair1
%load eigenpair2
% load eigenpair_NLS
% load eigenpair_NLS_v2 % Proof for Sec 5.1 & Sec 5.2
% load eigenpair_NLS_pt2 % Proof for Sec 5.3
% load eigenpair_CGL_pi_4

% CGL: pt1
% load ../Manifolds/eigenpair_CGL_0_pt1.mat
% filename = '../Manifolds/eigenpair_CGL_0_pt1_proof';
load ../Manifolds/eigenpair_CGL_pi_4_pt1.mat
filename = '../Manifolds/eigenpair_CGL_pi_4_pt1_proof';

% CGL: pt2
% load ../Manifolds/eigenpair_CGL_pi_4_pt2.mat
% filename = 'eigenpair_CGL_pi_4_pt2_proof';
% load ../Manifolds/eigenpair_CGL_0_pt2.mat
% filename = 'eigenpair_CGL_0_pt2_proof';

v = clean_p(v); x = clean_p(x);

nu = 1;

[I,success] = rad_poly_eigenpairs(x,v,theta,nu);

disp('   ')
disp(['I = [',num2str(I(1)),',',num2str(I(2)),']'])

m = (length(x)-1)/2; a = x(2:m+1); plot_periodic_complex(a)

lambda = x(1);
b = x(m+2:2*m+1);
r0 = I(1);

save(filename,'a', 'b', 'lambda', 'theta', 'r0', 'nu')
function [err,W_at_endpoint,W_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,d_all,a,h,N,n,theta,core,rigorous)
arguments
  eps_all; d_all; a; h; N; n; theta; core;
  rigorous = true;
end

%% start to solve variational problem
  [C,r_minus] = solve_variational_equation(a,theta,h,N,n,core,rigorous);
  % r_minus
  if rigorous>0
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(intval(C0(1,:,:))),2),2*core+1,1);
    M_phi = max(2*(sum(abs(intval(C)))'+r_minus)-C0k);
  else
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(C0(1,:,:)),2),2*core+1,1);
    M_phi = max(2*(sum(abs(C))'+r_minus)-C0k);
  end
  
  [C_backward,r_minus_backward] = solve_variational_equation_backward(a,theta,h,N,n,core,rigorous);
  if rigorous>0
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(intval(C0_backward(1,:,:))),2),2*core+1,1);
    M_psi = max(2*(sum(abs(intval(C_backward)))'+r_minus_backward)-C0k_backward);
  else
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(C0_backward(1,:,:)),2),2*core+1,1);
    M_psi = max(2*(sum(abs(C_backward))'+r_minus_backward)-C0k_backward);
  end
  
%   M0 = M_phi*M_psi;
% % A technique to estimate a uniform bound using interval arithmetic
  if rigorous
    Wm = getting_M0(C0,C0_backward,r_minus,r_minus_backward,h,n,core);%sup(Wm)
    Wmt = getting_M_t1s(C0,C0_backward,r_minus,r_minus_backward,h,n,core);
    Wm = min(M_phi*M_psi,Wm);
    Wmt = min(M_phi*M_psi,Wmt);
  else
    Wm = sup(getting_M0(C0,C0_backward,r_minus,r_minus_backward,h,n,core));
    Wm = min(M_phi*M_psi,Wm);
    Wmt = sup(getting_M_t1s(C0,C0_backward,r_minus,r_minus_backward,h,n,core));
    Wmt = min(M_phi*M_psi,Wmt);
  end
  
  C0(2:end,:,:) = 2*C0(2:end,:,:);
  if rigorous>0
    phi_at_endpoint = norm(reshape(sum(intval(C0),1),2*core+1,2*core+1),1);
  else
    phi_at_endpoint = norm(reshape(sum(C0,1),2*core+1,2*core+1),1);
  end
  
  %% start to estimate the evolution operator
  if rigorous>0
    ipi = intval('pi'); itheta = theta;
    mu_m = (core+1)^2*(2*ipi)^2*cos(itheta);
    ba_X = a_norm(intval(a));
    a_infty = a; a_infty(:,N+1) = 0;
    ba_inf_X = a_norm(intval(a_infty));
    ba_inf_X_dual = ba_inf_X;
  else
    theta = mid(theta);
    mu_m = (core+1)^2*(2*pi)^2*cos(theta);
    ba_X = a_norm(a);
    a_infty = a; a_infty(:,N+1) = 0;
    ba_inf_X = a_norm(a_infty);
    ba_inf_X_dual = ba_inf_X;
%     ba_inf_X_dual = a_norm(a_infty);
  end
  
  if mu_m>2*ba_X
    W_infinity = (1-exp(-(mu_m-2*ba_X)*h))/(mu_m-2*ba_X);
    W_bar_infty = (h-W_infinity)/(mu_m-2*ba_X);
    kappa = 1 - 4*Wm*ba_inf_X^2*W_bar_infty;
  else
    W_infinity = (exp((2*ba_X-mu_m)*h)-1)/(2*ba_X-mu_m);
    W_bar_infty = (W_infinity-h)/(2*ba_X-mu_m);
    kappa = 1 - 4*Wm*ba_inf_X^2*W_bar_infty;
  end
  W_infinite_sup = max(1,exp((2*ba_X-mu_m)*h));
%   M_infinite_at_endpoint = exp(-(mu_m-2*ba_X)*h);
  
  if kappa>0
    U_matrix = [Wm, 2*Wm*W_infinity*ba_inf_X_dual;...
      2*Wm*W_infinity*ba_inf_X, W_infinite_sup]/kappa;
    U_matrix_at_endpoint1 = [phi_at_endpoint+2*Wmt*ba_inf_X_dual*h*U_matrix(2,1),...
      2*Wmt*ba_inf_X_dual*h*U_matrix(2,2);...
      2*W_infinity*ba_inf_X*U_matrix(1,1),...
      exp((2*ba_X-mu_m)*h)+2*W_infinity*ba_inf_X*U_matrix(1,2)];
    U_matrix_J = [Wmt*(1+2*ba_inf_X_dual*h*U_matrix(2,1)),...
      2*Wmt*ba_inf_X_dual*h*U_matrix(2,2);...
      2*W_infinity*ba_inf_X*U_matrix(1,1),...
      W_infinite_sup+2*W_infinity*ba_inf_X*U_matrix(1,2)];
%     U_matrix = [Wm/kappa, 2*Wm*W_infinity*ba_inf_X/kappa;...
%       2*Wm*W_infinity*ba_inf_X/kappa, W_infinite_sup+4*Wm*W_infinity^2*ba_inf_X^2/kappa];
%     M = sup(norm(U_matrix,1));
%     U_matrix_at_endpoint1 = [phi_at_endpoint*(1+2*M_psi*ba_inf_X*h*U_matrix(2,1)),...
%       2*phi_at_endpoint*M_psi*ba_inf_X*h*U_matrix(2,2);...
%       2*W_infinity*ba_inf_X*U_matrix(1,1),...
%       exp((2*ba_X-mu_m)*h)+2*W_infinity*ba_inf_X*U_matrix(1,2)];
%     U_matrix_at_endpoint2 = [phi_at_endpoint*M_psi*(1+2*ba_inf_X*h*U_matrix(2,1)),...
%       2*phi_at_endpoint*M_psi*ba_inf_X*h*U_matrix(2,2);...
%       2*W_infinity*ba_inf_X*U_matrix(1,1),...
%       W_infinite_sup+2*W_infinity*ba_inf_X*U_matrix(1,2)];
%     M_at_endpoint = min(M,sup(norm(U_matrix_at_endpoint1,1)));
%     Ms = min(M,sup(norm(U_matrix_at_endpoint2,1)));
    W_h = sup(norm(U_matrix,1));
    F = @(x) W_h*(eps_all+h*(2*x.^2+d_all))-x;
    W_at_endpoint = min(W_h,sup(norm(U_matrix_at_endpoint1,1)));
    W_J = min(W_h,sup(norm(U_matrix_J,1)));
  else
    disp('Linearized problem is not solved (kappa<0)')
    W_at_endpoint = NaN; W_h = NaN; W_J = NaN; err = NaN;
    return
  end
  
  
  %% verify the contraction mapping: require INTLAB for using "verifynlssall"
%   F = @(x) M*(eps_all+stepsize*(2*x.^2+d_all))-x;
  
  [xx,xx_cand,data] = verifynlssall(F,infsup(0,1));
  if ~any(xx>0)
    disp('contraction mapping is not verified!')
    err = NaN;
    return
  end
  while(1)
    if isempty(xx_cand)
      err = sup(xx(:,all(F(sup(xx))<0,1)));
      break
    else
      [xx,xx_cand,data] = verifynlssall(data);
    end
  end
 
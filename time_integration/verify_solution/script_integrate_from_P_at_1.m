%% preliminary
clear
addpath('../verify_defect/')
addpath('../variational_problem/')

% load ../unstable_pt.mat
load ../end_points_proof1_CGL.mat

N = (size(P_at_1,1)-1); % # of Fourier projection
% N = 14; % # of Fourier projection
n = 13; % # of Chebyshev coefficients
stepsize = 0.08/2^5; % length of time step
tspan = [0,stepsize];
theta = intval('pi')/4;
% theta = pi/4;
rigorous = 1;

core = 2; % core = 2 is best?
a0 = [flipud(P_at_1(2:end));P_at_1].';
% a0 = zeros(1,2*N+1); % Initial sequence
% a0(N+1)=50; a0(N)=-25; a0(N+2)=-25;

if rigorous>0
  y_local = intval(zeros(1,10));
else
  y_local = zeros(1,10);
end

y = []; % Data container

%% getting approximate solution and residual bounds
% timestep = 1;
eps_all = error_at_1;
% timestep =
%     63
for timestep = 1:1
  [a, d_N, d_infty] = getting_the_solution_timestepping(N,n,tspan,a0,theta,rigorous);% Output is one-sided Chebyshev!
  disp(['delta_N = ',num2str(d_N)])
  disp(['delta_tail = ',num2str(d_infty)])
  h = stepsize;
  
%   plot_solution(a,tspan,1), hold on, pause(0.01)
  
%   if timestep==1
    % Initial error
    if rigorous>0
      eps_all = sup(eps_all + compute_eps(intval(a),intval(a0)));
      d_all = sup(d_N+d_infty);
    else
      eps_all = eps_all + compute_eps(a,a0);
      d_all = d_N+d_infty;
    end
%   end
    
%%%%%%%%%
  [err,M_at_endpoint,Ms,M0,M,a_X,kappa] = verify_local_existence(eps_all,d_all,a,h,N,n,theta,core);
  if any(isnan(err))
    break
  end
%%%%%%%%%%
  
  err_at_endpoint = intval(M_at_endpoint)*eps_all+intval(Ms)*stepsize*(2*intval(err)^2+intval(d_all));
  err_at_endpoint = min(err,sup(err_at_endpoint));
  
  %% verify global existence by Jonathan's method
  if rigorous > 0
    ipi = intval('pi'); itheta = theta;
    a0 = sum(intval(a))+sum(intval(a(2:end,:)));% initial sequence of next step
    if norm(a0,1) < err_at_endpoint
      error('global existence is failed...')
    end
    rs = norm(a0(1:N),1) + norm(a0(N+2:end),1) + err_at_endpoint;
    rc = abs(a0(N+1))+err_at_endpoint + 0.02*rs;
    ac = a0(N+1) + err_at_endpoint*infsup(-1,1);
%     theta = angle*ipi/180;
    mu = (2*ipi)^2*cos(theta);
    success_GE = verify_GE(rc,rs,mu,itheta,ac);
    a0 = mid(a0);
  else
    a0 = sum(a)+sum(a(2:end,:));% initial sequence of next step
    rs = norm(a0(1:N),1) + norm(a0(N+2:end),1) + err_at_endpoint;
    rc = abs(a0(N+1))+err_at_endpoint + 0.02*rs;
    ac = a0(N+1);
%     theta = angle/180*pi;
    mu = (2*pi)^2*cos(theta);
    success_GE = verify_GE(rc,rs,mu,theta,ac);
  end
  
  
  %% Data
  y_local(1) = tspan(1);
  y_local(2) = tspan(2);
  y_local(3) = M0;
  y_local(4) = M_at_endpoint;
  y_local(5) = M;
  y_local(6) = a_X;
  y_local(7) = kappa;
  y_local(8) = eps_all;
  y_local(9) = d_all;
  y_local(10) = err;
  % y(10) = sup(min(xx));
  y = [y;y_local];
  
  if success_GE>0
    break
  end
  
  %% Update the initial error and time interval
  eps_all = err_at_endpoint;
  tspan = tspan + stepsize;% next time step
  
end

success_GE
% printresult_timestepping
% save('data_GEv1.mat','y','t')
% save data_CGL_from_P_at_1.mat

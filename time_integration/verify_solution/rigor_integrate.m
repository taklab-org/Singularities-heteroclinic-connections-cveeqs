function [success,x,z,y] = rigor_integrate(n,P_at_1,error_at_1,theta,num_integration)
success = 0;
% P_at_1 should be one-sided Chebyshev
N = (size(P_at_1,1)-1); % # of Fourier projection
stepsize = 0.08/2^6; % length of time step
tspan = [0,stepsize];


core = 0; % core = 2 is best?
a0 = mid([flipud(P_at_1(2:end));P_at_1].');
% error_at_1 = 0;
err_at_endpoint_old = error_at_1;
% eps_all = error_at_1;

x = zeros(n*num_integration,1);
z = zeros(n*num_integration,1);
% x = [];
% z = [];

y_local = intval(zeros(1,10));
y = []; % Data container

%% getting approximate solution and residual bounds
% timestep = 1;

% timestep =
%     63

% initial_step = true;

% for timestep = 1:num_integration
timestep = 1;
while timestep <= num_integration
  disp(['timestep: ',num2str(timestep)])
  [a, d_N, d_infty] = getting_the_solution_timestepping(N,n,tspan,a0,theta,1);% Output is one-sided Chebyshev!
%   disp(['delta_N = ',num2str(d_N)])
%   disp(['delta_tail = ',num2str(d_infty)])
  h = stepsize;
  
  eps_all = sup(err_at_endpoint_old + compute_eps(intval(a),intval(a0)));
  d_all = sup(d_N+d_infty);
  
  %%%%%%%%%
  [err,M_at_endpoint,Ms,M0,M,a_X,kappa] = verify_local_existence(eps_all,d_all,a,h,N,n,theta,core);
  if any(isnan(err))
    disp('local existence is failed!!')
    break
  end
  %%%%%%%%%%
  
  err_at_endpoint = intval(M_at_endpoint)*eps_all+intval(Ms)*stepsize*(2*intval(err)^2+intval(d_all));
  err_at_endpoint = min(err,sup(err_at_endpoint));
  disp(['error at endpoint = ',num2str(err_at_endpoint)])
  
%   if timestep==430
%       pause
%   end
  
  % adjust stepsize corresponding to increase ratio of err_at_endpoint
  if any(err_at_endpoint./err_at_endpoint_old>1.05) && timestep~=1
    stepsize = stepsize/(max(err_at_endpoint./err_at_endpoint_old))^2;
    tspan(2) = tspan(1) + stepsize;% next time step
    core = max(0, core-1);
    disp('adjust timestep (smaller)')
    continue
  elseif all(err_at_endpoint./err_at_endpoint_old <= 1.01) && (stepsize < 2.5e-3)% max stepsize
    stepsize = stepsize*1.1;
    tspan(2) = tspan(1) + stepsize;% next time step
    core = min(2, core+1);
    disp('adjust timestep (larger)')
    continue
  end
  
  %% verify global existence by Jonathan's method
  ipi = intval('pi'); itheta = theta;
  a0 = sum(intval(a))+sum(intval(a(2:end,:)));% initial sequence of next step
  if norm(a0,1) < err_at_endpoint
    error('global existence is failed...')
  end
  rs = norm(a0(1:N),1) + norm(a0(N+2:end),1) + err_at_endpoint;
  rc = abs(a0(N+1))+err_at_endpoint + 0.02*rs;
  ac = a0(N+1) + err_at_endpoint*infsup(-1,1);
  mu = (2*ipi)^2*cos(itheta);
  success = verify_GE(rc,rs,mu,itheta,ac);
  a0 = mid(a0);
  
%   initial_step = false;

  %   %% Data
  y_local(1) = tspan(1);
  y_local(2) = tspan(2);
  y_local(3) = M0;
  y_local(4) = M_at_endpoint;
  y_local(5) = M;
  y_local(6) = a_X;
  y_local(7) = kappa;
  y_local(8) = err_at_endpoint;
  y_local(9) = d_all;
  y_local(10) = err;
  % y(10) = sup(min(xx));
  y = [y;y_local];

  uval=ifft(ifftshift(chebcoeffs2chebvals([a(1,:);2*a(2:end,:)]),2).');
  % plot(chebpts(n,tspan),max(abs(y)*size(y,1)))
  z((timestep-1)*n+1:timestep*n) = max(abs(uval)*size(uval,1));
%   if angle_count == 1
  x((timestep-1)*n+1:timestep*n) = chebpts(n,tspan);
%   end
  
  if success>0
    break
  end
  
  %% Update the initial error and time interval
  err_at_endpoint_old = err_at_endpoint;
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + stepsize;% next time step
  timestep = timestep + 1;
end

zero_index = find(x(2:end)==0,1,'first');
x(zero_index:end) = []; z(zero_index:end) = [];

% plot(mid(y(:,2)),mid(y(:,8)),'Linewidth',1.6)

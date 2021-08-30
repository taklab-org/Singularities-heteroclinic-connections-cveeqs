%%
clear
% clf
% addpath('../../toolbox/chebfun-master/')
% addpath('../verify_solution/')
% addpath('../verify_defect/')
addpath('../../Manifolds/')

n = 30; % # of Chebyshev coefficients
stepsize = 0.08/2^4; % length of time step
rigorous=0;
gamma = exp(1i*pi/4);
% gamma = 1i;
% gamma = 1;
num_integration = 50;
% sol_data = zeros(n*num_integration,1,n*num_integration);
% n_pad = 10;

% load unstable_pt.mat
for ii=328
  clf
  theta = ii/180*pi;
  [P_at_minus_1,P_at_1] = produce_unstable_data_from_angle_for_picture(theta,10);
%   if (ii > 200) && (ii < 260)
%       n_pad = 20;
%   else
      n_pad = 0;
%   end
%   a0 = mid([flipud(P_at_1(2:end));P_at_1].');
  a0 = mid([zeros(n_pad,1);flipud(P_at_1(2:end));P_at_1;zeros(n_pad,1)].');

  
  ii
  
  N = (size(a0,1)-1)/2; % # of Fourier coefficients
    
%   figure
  plot_profile(a0,['r','b']),hold on

  %% getting approximate solution
  % timestep = 1;
  chebfunpref.setDefaults('factory');
  opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
  tspan = [0,stepsize];
  for timestep = 1:num_integration
    chebfunpref.setDefaults('fixedLength',n);
    %
    u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,a0,opts);
    % u_cheb = chebfun.ode113(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
    % u_cheb = chebfun.ode15s(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
    ta = chebcoeffs(u_cheb);
%     y=ifft(ifftshift(chebcoeffs2chebvals(ta),2).');
%     plot(chebpts(n,tspan),max(abs(y)*size(y,1)))

    ta = [ta(1,:);0.5*ta(2:end,:)];
    ta_fourier = sum(ta);
    if abs(ta_fourier(end))>0.99
      break
    end
%     plot_solution(ta,tspan,1),hold on, pause(0.01)
    if mod(timestep,5)==4, plot_profile(a0,'k'),hold on, pause(0.01),end
    a0 = u_cheb(end);% initial sequence of next step
    tspan(1) = tspan(2);
    tspan(2) = tspan(2) + stepsize;% next time step
  end
  chebfunpref.setDefaults('factory');
  sgtitle(['CGL, angle: ',num2str(ii)])
%   SaveFig(gcf,['figs/CGL_angle',num2str(ii)])
  
end
function dy = ode_func(y,gamma)

N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;

end


%%
clear
clf
addpath('../toolbox/chebfun-master/')
addpath('verify_solution/')
addpath('verify_defect/')

load unstable_pt2.mat

N = (size(a_unstable2,1)-1)/2; % # of Fourier coefficients
n = 13; % # of Chebyshev coefficients
stepsize = 0.08/2^4; % length of time step
tspan = [0,stepsize];
rigorous=0;
% theta = 45/180*pi;
gamma = exp(1i*theta);


a0 = a_unstable2.';
% a0 = zeros(1,2*N+1); % Initial sequence
% a0(N+1)=50; a0(N)=-25; a0(N+2)=-25;


% plot_profile(a0,'r'),hold on

%% getting approximate solution
% timestep = 1;
chebfunpref.setDefaults('factory');
opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
for timestep = 1:20
chebfunpref.setDefaults('fixedLength',n);
%
% opts = odeset('abstol',1e-18,'reltol',1e-18);
u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,a0,opts);
% u_cheb = chebfun.ode113(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% u_cheb = chebfun.ode15s(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
ta = chebcoeffs(u_cheb);
ta = [ta(1,:);0.5*ta(2:end,:)];
% plot_solution(ta,tspan,1),hold on, pause(0.01)
% if mod(timestep,5)>0, plot_profile(a0,['b','r']),hold on, pause(0.1),end
a0 = u_cheb(end);% initial sequence of next step
tspan = tspan + stepsize;% next time step
% plot(tspan(2),abs(real(a0(N+1))),'x'), hold on, pause(0.1)
plot_solution(ta,tspan,1), hold on, pause(0.1)
end

chebfunpref.setDefaults('factory');
% plot_profile(a_unstable2.','k'),hold on

function dy = ode_func(y,gamma)

N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;

end


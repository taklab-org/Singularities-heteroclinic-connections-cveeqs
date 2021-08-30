%%
clear
% clf
% addpath('../../toolbox/chebfun-master/')
addpath('verify_solution/')
addpath('verify_defect/')

addpath('../Manifolds/')


% load unstable_pt.mat
for ii=180
theta = ii/180*pi;
[P_at_minus_1,P_at_1] = produce_unstable_data_from_angle(theta,10);
% a0 = mid([flipud(P_at_minus_1(2:end));P_at_minus_1].');
a0 = mid([flipud(P_at_1(2:end));P_at_1].');

ii

N = (size(a0,1)-1)/2; % # of Fourier coefficients
n = 30; % # of Chebyshev coefficients
stepsize = 0.08/2^4; % length of time step
tspan = [0,stepsize];
rigorous=0;
% gamma = exp(1i*pi/4);
gamma = 1i;




% a0 = a_unstable.';
% a0 = zeros(1,2*N+1); % Initial sequence
% a0(N+1)=50; a0(N)=-25; a0(N+2)=-25;

figure
plot_profile(a0,['r','b']),hold on

%% getting approximate solution
% timestep = 1;
chebfunpref.setDefaults('factory');
opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
for timestep = 1:100
chebfunpref.setDefaults('fixedLength',n);
%
% opts = odeset('abstol',1e-18,'reltol',1e-18);
u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,a0,opts);
% u_cheb = chebfun.ode113(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% u_cheb = chebfun.ode15s(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
ta = chebcoeffs(u_cheb);
ta = [ta(1,:);0.5*ta(2:end,:)];
% plot_solution(ta,tspan,3),hold on, pause(0.01)
if mod(timestep,5)==4 plot_profile(a0,'k'),hold on, pause(0.01),end
a0 = u_cheb(end);% initial sequence of next step
tspan = tspan + stepsize;% next time step
% [abs(a0(N+1)),norm([a0(1:N);a0(N+2:end)],1)]
% if abs(a0(N+1))<0.176
%   return
% end
% plot(tspan(2),abs(real(a0(N+1))),'x'), hold on
% pause
end

chebfunpref.setDefaults('factory');


% phi=a0(N+1);
% % phi = -real(phi)+1i*imag(phi);
% L = chebop([0,5]);
% L.op = @(t,y) diff(y)-2*phi/(1-gamma*phi*t)*y;
% L.lbc = phi;
% y = L\0;
% 
% % plot(abs(y))
% plot(y)
end
function dy = ode_func(y,gamma)

N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;

end


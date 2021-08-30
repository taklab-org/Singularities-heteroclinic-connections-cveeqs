%%
clear
% clf
addpath('verify_solution/')
addpath('verify_defect/')
addpath('variational_problem_NLS/')


num_integration = 1e4;

% angle at eigen vector:
% angle = 0;
% angle = [0:179,181:359];
angle = [0:18:359,359];
% angle = 0:359;
% angle = angle(mod(angle,4)==0);
% x = zeros(n*num_integration,1);
angle_eigvec = angle/180*pi;
% z = zeros(n*num_integration,length(angle));
angle_count = 1;

% angle of complex plain of time:
% theta = intval('pi')/4;
% theta = 0;

for ii = angle
  % Produce unstable manifolds
  addpath('../Manifolds/')
  [P_at_1,error_at_1] = produce_unstable_data_from_angle_NLS(ii/180*pi,10);
  rmpath('../Manifolds/')
  n_pad = 10;
  a0 = mid([zeros(n_pad,1);flipud(P_at_1(2:end));P_at_1;zeros(n_pad,1)].');
    
  N = (size(a0,1)-1)/2; % # of Fourier coefficients
  
  disp(['angle at eig v.: ', num2str(ii)])
    
    
  %% rigorous integration
  [success,x,z,data] = rigor_integrate_NLS(P_at_1,error_at_1,num_integration);
  y = angle_eigvec(angle_count)*ones(size(x));
  surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], 'EdgeColor','flat', 'FaceColor','none','Linewidth',2),pause(0.01),hold on
  if success==0
      plot(x(end),y(end),'rx','Markersize',16,'Linewidth',2)
  end
  save(['data_NLS_angle',num2str(ii)],'data','x','y','z','success')
  angle_count = angle_count + 1;
end

view(3)
% [X,Y] = meshgrid(x,y);
% surf(X,Y,z','EdgeColor','none')
colorbar
xlabel('$t$','interpreter','latex')
ylabel('Angle at eigen vec.','interpreter', 'latex')
% zlabel('$\|\bar{u}\|_{\infty}$','interpreter', 'latex')
% title('Solution distribution v.s. angle')
view([0, 90])
set(gca,'ColorScale','log')
caxis([1 8e2])
xlim([0,0.2])
ylim([0,2*pi])
% % SaveFig(gcf,'figs/CGL_0_pt2_solution_distribution')



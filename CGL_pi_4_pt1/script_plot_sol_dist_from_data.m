% Plot solution distribution given by rigor_sol_dist_CGL.m
% dataFiles = dir('data_CGL_0_angle*');
dataFiles = dir('data_CGL_pi_4_angle*');
numfiles = length(dataFiles);

for k = 1:numfiles
  clear x y z success data
  load(dataFiles(k).name) 
  surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], 'EdgeColor','flat', 'FaceColor','none','Linewidth',2);hold on
  if success==0
    disp(dataFiles(k).name)
      plot3(x(end),y(end),10000,'rx','Markersize',8,'Linewidth',2);
  end
end

view(3)

colorbar
% xlabel('$t$','interpreter','latex')
% ylabel('Angle at eigen vec.','interpreter', 'latex')
% zlabel('$\|\bar{u}\|_{\infty}$','interpreter', 'latex')
% title('Solution distribution v.s. angle')
view([0, 90])
set(gca,'ColorScale','log')
caxis([1 8e2])
xlim([0,0.2])
ylim([0,2*pi])
set(gca,'FontSize',20)
xlabel('$t$','interpreter','latex','FontSize', 30)
ylabel('Angle at eigen vector $\psi_k$','interpreter', 'latex','FontSize', 30)
hcb=colorbar;clabel = get(hcb,'Label');titleString = '$\|u(t)\|_{\infty}$';set(clabel ,'String',titleString,'interpreter', 'latex');
yticks([0 pi/2 pi 3*pi/2 2*pi]),yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
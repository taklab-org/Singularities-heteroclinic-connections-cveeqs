% Plot solution distribution given by rigor_sol_dist_CGL.m
% dataFiles = dir('data_CGL_0_angle*');
dataFiles = dir('data_CGL_pi_4_angle*');
numfiles = length(dataFiles);

for k = 1:numfiles
  clear x y z success data
  load(dataFiles(k).name) 
  surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], 'EdgeColor','flat', 'FaceColor','none','Linewidth',2);hold on
  if success==0
      plot(x(end),y(end),'rx','Markersize',16,'Linewidth',2);
  end
end

view(3)
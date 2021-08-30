clc
clear
% close all

% tot_eq contains five steady states 
% computed previously with Newton's method
load tot_eq 

N = 14;
theta = 0;

a1 = tot_a(:,5);% 5 gives pt1, 3 gives pt2
% a1 = [a1;zeros(20,1)];
a1=a1(1:N+1);
% a2 = conj(a1); a2 = newton(a2,theta);

DF1 = DF_steady_states(a1,theta);
% plot_periodic_complex(a1)
[V,E] = eig(DF1);
eigenvalues = diag(E);
% ind = (real(eigenvalues)>0);
[~,ind] = max(real(eigenvalues));
unstable_eigs1 = eigenvalues(ind)
unstable_eigvectors1 = V(:,ind);

a = a1; b = unstable_eigvectors1(:,1); lambda = unstable_eigs1(1);
x = [lambda;a;b]; v = b; 
% save eigenpair2_NLS x v theta
% save ../Manifolds/eigenpair_CGL_pi_4_pt1 x v theta

% figure
% DF2 = DF_steady_states(a2,theta);
% plot_periodic_complex(a2)
% [V,E] = eig(DF2);
% eigenvalues = diag(E);
% [~,ind] = max(real(eigenvalues));
% unstable_eigs2 = eigenvalues(ind)
% unstable_eigvectors2 = V(:,ind);
% 
% 
% a = a2; b = unstable_eigvectors2(:,1); lambda = unstable_eigs2(1);
% x = [lambda;a;b]; v = b; 
% % save eigenpair_NLS_pt2_conj x v theta
% Plot figure
plot(eigenvalues(real(eigenvalues)<=0),'x','MarkerSize',18,'linewidth',3)
hold on
plot(eigenvalues(real(eigenvalues)>0),'*','MarkerSize',18,'linewidth',3)
set(gca,'FontSize',20)
xlabel('$Re\,(\tilde{\lambda})$','interpreter','latex','FontSize', 30)
ylabel('$Im\,(\tilde{\lambda})$','interpreter', 'latex','FontSize', 30)
% plot([-8000,1000],[0,0],'k--','linewidth',2)
% plot([0,0],[-40,40],'k--','linewidth',2)
hold off

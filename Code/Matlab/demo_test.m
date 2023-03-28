addpath("./func");

f = @(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
g = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), -200*(x(1)^2-x(2))]';
x0 = [-1000 100]';

[x_ls,df_ls,x] = show_bfgs(f,g,x0);
% figure(); hold on;
semilogy(df_ls);
xlabel('Index','FontSize',16,'interpreter','latex');
ylabel('$\nabla f$','FontSize',16,'interpreter','latex');
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));
% set()
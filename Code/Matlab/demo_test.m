addpath("./func");

f = @(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
g = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), -200*(x(1)^2-x(2))]';
x0 = [-1000 100000]';

[x,cal_f,cal_g] = bfgs(f,g,x0);

% [x_ls,df_ls,x] = show_bfgs(f,g,x0);
% % figure(); hold on;
% semilogy(df_ls);
% xlabel('Index','FontSize',16,'interpreter','latex');
% ylabel('$||\nabla f||$','FontSize',16,'interpreter','latex');
% grid on; box on;
% set(gca,'Color','none');
% set(gca,'LooseInset',get(gca,'TightInset'));
% set()


% x = -1000:0.5:1000;
% y = -1000:0.5:1000;
% [X,Y] = meshgrid(x,y);
% Z = 100*(X.^2-Y).^2+(X-1).^2;
% contour(X,Y,Z,20)
% 
% a=1;
% b=100;
% xa=-3:0.1:3; ya=-3:0.1:3;
% [x,y]=meshgrid(xa,ya);
% z=(a-x).^2+b*(y-x.^2).^2;
% mesh(x,y,z)
% f = @(x,y) 100*(x.^2-y).^2+(x-1).^2;
% fcontour(f,[-abs(x0(1)) abs(x0(1)) -abs(x0(2)) abs(x0(2))],20)

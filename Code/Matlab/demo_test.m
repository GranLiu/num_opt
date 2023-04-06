% -------------------------------------------------------------------------
% A test demo for BFGS, using Rosenbrock optimization problem.
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch6.1]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch5.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
clc;
clear;
close all;

addpath("./func");

%% Rosenbrock function
f = @(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
g = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1), -200*(x(1)^2-x(2))]';
x0 = [10000 10000]';

%% use show_* to extract gradient and point at each itr
% [x,cal_f,cal_g,itx] = bfgs(f,g,x0);
[x_ls,cal_f,cal_g,itx,df_ls,x1] = show_bfgs(f,g,x0);
disp(['Count of f(`): ',num2str(cal_f),' Count of g(`): ',num2str(cal_g)])

%% plot gradient
figure()
semilogy(df_ls);
xlabel('Index','FontSize',16,'interpreter','latex');
ylabel('$||\nabla f||$','FontSize',16,'interpreter','latex');
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));

%% plot function value
f_val = zeros(1,length(x_ls));
for idx = 1:length(x_ls)
    f_val(idx) = f(x_ls(:,idx));
end
figure()
% plot(1:length(x_ls),f_val);
semilogy(f_val);
xlabel('Index','FontSize',16,'interpreter','latex');
ylabel('$f$','FontSize',16,'interpreter','latex');
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));

%% plot contours and surface
x = -max(x_ls,[],'all'):0.1:max(x_ls,[],'all');
y = -max(x_ls,[],'all'):0.1:max(x_ls,[],'all');
figure()
[X,Y] = meshgrid(x,y);
Z = 100*(X.^2-Y).^2+(X-1).^2;
[C,h] = contour(X,Y,Z,100); hold on;
clabel(C,h);
plot(x_ls(1,:),x_ls(2,:),'k-o','linewidth',1.4);
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));

figure()
xa=-max(x_ls,[],'all'):0.1:max(x_ls,[],'all');
ya=-max(x_ls,[],'all'):0.1:max(x_ls,[],'all');
[x,y]=meshgrid(xa,ya);
z=(1-x).^2+100*(y-x.^2).^2;
mesh(x,y,z);
hold on;
plot3(x_ls(1,:),x_ls(2,:),f_val,'k-o','linewidth',1.4);
xlabel('$x_1$','interpreter','latex','fontsize',16);
ylabel('$x_2$','interpreter','latex','fontsize',16);
zlabel('$f(\mathbf{x})$','interpreter','latex','fontsize',16);
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));

figure()
s = surfc(x,y,z);
zlim([-0.5*max(z,[],'all'),max(z,[],'all')]); hold on;
plot3(x_ls(1,:),x_ls(2,:),-0.5*max(z,[],'all')*ones(1,length(x_ls)),'g-*','linewidth',1.6);
plot3(x_ls(1,:),x_ls(2,:),f_val,'g-o','linewidth',1.6);
s(2).LineWidth = 2;
s(2).LevelStep = 100;
grid on; box on;
set(gca,'Color','none');
set(gca,'LooseInset',get(gca,'TightInset'));
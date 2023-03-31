function [a,cal_f,cal_g] = ls_aj(f,g,x0,p,a,c1)
% -------------------------------------------------------------------------
% Inexact line search to calculate step length a -- Armijo condition.
% Use quadratic interpolation to calculate a at each iteration.
% 
% Input
%   f: objective function
%   g: gradient of objective function
%   x0: start point
%   p: search direction
%	a: intial step length
%	c1: slope for Armijo condition [DEF: 1e-4 \in (0,1)]
% Output
%   a: step length satisfying Wolfe condition
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch3.5]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch2.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
if nargin<4 || nargin>6
    error("There should be 4~6 inputs in this func.");
end
if nargin == 4
    a = 1;
    c1 = 1e-4;
elseif nargin == 5
    c1 = 1e-4;
end

phi = @(a) f(x0+a*p);

phi_0 = feval(phi,0);
g_val = feval(g,x0);
g_phi_0 = g_val.'*p;

% calculation count of f_val & g_val
cal_f = 1;
cal_g = 1;

a0 = a;
if phi(a0) <= phi_0+c1*a0*g_phi_0
    cal_f = cal_f + 1;
    return;
end
a1 = -g_phi_0*a0^2 / ( 2* (phi(a0)-phi_0-g_phi_0*a0) );
cal_f = cal_f + 1;

% if a1/a0<1e-2
%     a1 = 0.1*a0;
% end

if phi(a1) <= phi_0+c1*a1*g_phi_0
    a = a1;
    cal_f = cal_f + 1;
    return;
end
tmp_a1 = phi(a1)-phi_0-g_phi_0*a1;
tmp_a0 = phi(a0)-phi_0-g_phi_0*a0;
while phi(a) > phi_0+c1*a*g_phi_0
    cal_f = cal_f + 1;
    zz = 1/(a0^2*a1^2*(a1-a0)) * [a0^2 -a1^2; -a0^3 a1^3] *...
        [tmp_a1; tmp_a0];
    a0 = a1;
    a1 = ( -zz(2) + sqrt(zz(2)^2-3*zz(1)*g_phi_0) ) / (3*zz(1));
    % avoid either too close or too smaller than its prodecessor, p58 in [1]
%     if a1/a0<0.1 || a1/a0>0.9
%         a1 = 0.5*a0;
%     end
    tmp_a0 = tmp_a1;
    tmp_a1 = phi(a1)-phi_0-g_phi_0*a1;
    cal_f = cal_f + 1;
    a = a1;
end
end
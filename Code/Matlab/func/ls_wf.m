function a = ls_wf(f,g,x0,p,a,c1,c2)
% -------------------------------------------------------------------------
% Inexact line search to calculate step length a -- Wolfe condition.
%
% Input
%   f: objective function
%   g: gradient of objective function
%   x0: start point
%   p: search direction
%	a: intial step length
%	c1: slope for Armijo condition [DEF: 1e-4 \in (0,1)]
%   c2: slope for curvature condition [DEF: 0.9 for Newton, 0.01 for CG \in (c1,1)]
% Output
%   a: step length satisfying Wolfe condition
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch3.5]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch2.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
if nargin<4 || nargin>7
    error("There should be 4~7 inputs in this func.");
end
if nargin == 4
    a = 1;
    c1 = 1e-4;
    c2 = 0.9;
elseif nargin == 5
    c1 = 1e-4;
    c2 = 0.9;
end
% Algorithm 3.5 in [1], page 60
a0 = 0;
a1 = a;
a_max = 3;
phi = @(a) f(x0+a*p);
g_phi = @(a) g(x0+a*p).'*p;
phi_0 = phi(0);
g_phi_0 = g_phi(0);
phi_a0 = phi(a0);
idx = 1;
while true
    phi_a1 = phi(a1);
    if phi_a1>phi_0+c1*a1*g_phi_0 || (phi_a1>=phi_a0 && idx>1)
        a = zoom(phi,g_phi,a0,a1,c1,c2);
        break;
    end
    phi_a0 = phi_a1;
    d_phi_a1 = g_phi(a1);
    if abs(d_phi_a1) <= -c2*g_phi_0
        a = a1;
        break;
    end
    if d_phi_a1>=0
        a = zoom(phi,g_phi,a1,a0,c1,c2);
        break;
    end
    a0 = a1;
    a1 = a0 + (a_max-a0) / 4;
    idx = idx +1;
end

end
function [a,cal_f,cal_g] = zoom(phi,g_phi,a_l,a_h,c1,c2)
% -------------------------------------------------------------------------
% Zoom the interval in calculatation of step length a -- Wolfe condition.
%
% Input
%   phi: eq.3.54 in [1]
%   g_phi: left-hand-side of eq.3.5, gradient of phi
%	a_l: step length lower bound
%   a_h: step length higher bound
%	c1: slope for Armijo condition
%   c2: slope for curvature condition
% Output
%   a: step length satisfying Wolfe condition
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch3.5]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch2.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
if nargin~=6
    error("6 params required in ZOOM!");
end

itx = 1;
itx_max = 20;

phi_0 = phi(0);
g_phi_0 = g_phi(0);
cal_f = 1;
cal_g = 1;
while itx<itx_max
    % use bisection to find a trial step
    a = 0.3*(a_l+a_h);
    phi_a = phi(a);
    cal_f = cal_f + 1;
    if phi_a > phi_0+c1*a*g_phi_0 || phi_a >= phi(a_l)
        cal_f = cal_f + 1;
        a_h = a;
    else
        d_phi_a = g_phi(a);
        cal_g = cal_g + 1;
        if abs(d_phi_a) <= -c2*g_phi_0
            break;
        end
        if d_phi_a*(a_h-a_l) >= 0
            a_h = a_l;
        end
        a_l = a;
    end
    itx = itx+1;
end

end
function [a,cal_f,cal_g] = ls_bt(f,g,x0,p,a,rou,c1)
% -------------------------------------------------------------------------
% Inexact line search -- Backtracking. [Armijo Condition]
% Well-suited for Newton method but is less appropriate for quasi-Newton/CG
%
% Input
%   f: objective function
%   g: gradient of objective function
%   x0: start point
%   p: search direction
%	a: intial step length [DEF: 1]
%   rou: backtracking ratio [DEF: 0.618]
%	c: slope for Armijo condition [DEF: 1e-4 \in (0,1)]
% Output
%   a: step length satisfying Armijo condition
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch3.1]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch2.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
idx = 0;
idx_max = 20;
if nargin<4 || nargin>7
    error("There should be 4~7 inputs in this func.");
end
if nargin == 4
    a = 1;
    rou = 0.618;
    c1 = 1e-4;
elseif nargin == 5
    rou = 0.618;
    c1 = 1e-4;
elseif nargin == 6
    c1 = 1e-4;
end
f_val = feval(f,x0);
g_val = feval(g,x0);
tmp = g_val.'*p;
% calculation count of f_val & g_val
cal_f = 1;
cal_g = 1;

while idx < idx_max
    cal_f = cal_f + 1;
    if feval(f,x0+a*p) <= f_val+c1*a*tmp
        break;
    end
    a = rou*a;
    idx = idx + 1;
end
end
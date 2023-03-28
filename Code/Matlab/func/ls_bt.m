function a = ls_bt(f,g,x0,p,a,rou,c)
% -------------------------------------------------------------------------
% Inexact line search to calculate step length a -- Backtracking.
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
if nargin<4 || nargin>7
    error("There should be 4~7 inputs in this func.");
end
if nargin == 4
    a = 1;
    rou = 0.618;
    c = 1e-4;
elseif nargin == 5
    rou = 0.618;
    c = 1e-4;
elseif nargin == 6
    c = 1e-4;
end
f_val = feval(f,x0);
g_val = feval(g,x0);
tmp = g_val.'*p;
while( feval(f,x0+a*p) > f_val+c*a*tmp )
    a = rou*a;
end
end
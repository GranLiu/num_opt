function [x_ls,df_ls,x1] = show_bfgs(f,g,x0,tor,h)
% -------------------------------------------------------------------------
% Quasi-Newton method -- BFGS
%
% Input
%   f: objective function
%   g: gradient of objective function
%	x0: starting point
%   tor: convergence tolerance [DEF: 1e-5]
%	h0: inverse Hessian approxiamtion [DEF: I]
% Output
%   x1: minimizer given by BFGS method
%
% Reference
% [1] "Numerical Optimization" -- Jorge Nocedal, Stephen J.Wright [Ch6.1]
% [2] "最优化方法及其 Matlab 程序设计" -- 马昌凤 [Ch5.2]
%
% Yongxi Liu, Xi'an Jiaotong University, 2023-03.
% -------------------------------------------------------------------------
if nargin<3 || nargin>5
    error("There should be 3~5 inputs in this func.");
end
len = length(x0);
if nargin == 3
    tor = 1e-5;
    h = eye(len);
elseif nargin == 4
    h = eye(len);
end
x1 = x0;
df = sum(abs(feval(g,x1)));
x_ls = [x0];
df_ls = [df];
while( df > tor )
    p = -h*feval(g,x1);
    x0 = x1;
%     a = ls_wf(f,g,x0,p);
    a = ls_aj(f,g,x0,p);
    a = ls_bt(f,g,x0,p);
    x1 = x0 + a*p;
    df = sum(abs(feval(g,x1)));
    x_ls = [x_ls x1];
    df_ls = [df_ls df];
    s = a*p;
    y = feval(g,x1) - feval(g,x0);
    rou = 1/(y.'*s);
    h = (eye(len) - rou*s*y.')*h*(eye(len)-rou*y*s.') + rou*(s*s.');
end

end
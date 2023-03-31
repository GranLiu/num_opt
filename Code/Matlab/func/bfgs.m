function [x1,cal_f,cal_g] = bfgs(f,g,x0,tor,h)
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
cal_f = 0;
cal_g = 0;
while( sum(abs(feval(g,x1))) > tor )
    % to avoid calculate gradient of fk again later (compromise spa->time)
    g_fk = feval(g,x1);
    p = -h*g_fk;
    x0 = x1;
%     [a,cal_f_tmp,cal_g_tmp] = ls_wf(f,g,x0,p);
    [a,cal_f_tmp,cal_g_tmp] = ls_aj(f,g,x0,p);
%     [a,cal_f_tmp,cal_g_tmp] = ls_bt(f,g,x0,p);
    
    % PROBLEM: if |s| is too short, y \aprx 0, leading to rou = INF
    s = a*p;
    x1 = x0 + s;
    y = feval(g,x1) - g_fk;
    rou = 1/(y.'*s);
    h = (eye(len) - rou*s*y.')*h*(eye(len)-rou*y*s.') + rou*(s*s.');
    
    % count of the f/g calculation
    cal_g = cal_g + cal_g_tmp + 2;
    cal_f = cal_f + cal_f_tmp;
end

end
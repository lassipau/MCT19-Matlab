function sol = LinSysSim(A,B,x0,ufun,tspan)
% function sol = LinSysSim(A,B,x0,ufun,tspan)
%
% Simulate the state of the differential equation x'(t)=Ax(t)+Bu(t)
% with initial state x(0)=x0, and u(t) = ufun(t) ('ufun' is a function
% handle) over the time interval 'tspan'. The returned variable 'sol' is
% the output of the Matlab's differential equation solver 'ode15s'.

odefun = @(t,x) A*x + B*ufun(t);

sol = ode15s(odefun,tspan,x0);
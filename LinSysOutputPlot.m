function [tt,yy] = LinSysOutputPlot(sol,C,D,ufun,N,axlim,LineW)
% function [tt,yy] = LinSysOutputPlot(sol,C,D,ufun,N)
%
% Plots the measured output of a linear system when 'sol' is the solution
% variable obtained from the ODE solver, C and D are parameters of the 
% system and 'ufun' is the function handle for the input function. Uses a 
% uniform grid with N points. 
% 'axlim' are the limits for the axes (input '[]' for default) and 'LineW' 
% is the line width.
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)


tt = linspace(sol.x(1),sol.x(end),N);
yy = C*deval(sol,tt)+D*ufun(tt);

if nargin <= 6
    LineW = 1;
end


plot(tt,yy,'Linewidth',LineW);

if nargin >5 && ~isempty(axlim)
    axis(axlim)
end
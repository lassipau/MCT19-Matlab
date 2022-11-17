function [tt,xx] = LinSysStatePlot(sol,N,axlim,LineW)
% function [tt,xx] = LinSysStatePlot(sol,N)
%
% Plots the state variables of a linear system when 'sol' is the solution
% variable obtained from the ODE solver. Uses a uniform grid with N points.
% 'axlim' are the limits for the axes (input '[]' for default) and 'LineW' 
% is the line width.
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

tt = linspace(sol.x(1),sol.x(end),N);
xx = deval(sol,tt);

if nargin <= 3
    LineW = 1;
end

plot(tt,xx,'Linewidth',LineW);

if nargin >2 & ~isempty(axlim)
    axis(axlim)
end
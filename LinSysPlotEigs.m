function LinSysPlotEigs(A,axlim)
% function PlotEigs(A,axlims)
%
% Plots the eigenvalues of A
% If 'axlim' is not given, limits determined from the spectrum. 
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

Aspec = eig(full(A));

if nargin == 1
  axlim = [min(real(Aspec)), max(real(Aspec)), min(imag(Aspec)), max(imag(Aspec))];
end
  

hold off
cla
hold on
plot(real(Aspec),imag(Aspec),'r.','Markersize',15)
% set the limits of the plot
axis(axlim)

% plot the axes
plot(axlim(1:2),[0 0],'k',[0 0],axlim(3:4),'k','Linewidth',1)

maxreal = num2str(max(real(Aspec)));

title(['Largest real part = $' maxreal '$' ],'Interpreter','Latex')
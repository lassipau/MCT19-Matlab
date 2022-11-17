% Example: Proportional-Integral control for the damped harmonic oscillator.
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

r = 1; k = 1; m = 1;


A = [0 1;-k/m -r/m];
B = [0;1/m];
C = [1 0];


% Construct the PI-controller
% Choose parameters K_P to stabilize A+B*K_P*C, 
% and the gain parameter eps>0

% "Root locus" can be used to investigate the effect of K_P on the
% eigenvalues of A+B*K_P*C (the two lines commented out). 
% Here we see that the real parts of the matrix do not depend on K_P, 
% so the stability margin does not improve with K_P. We can choose K_P=0.
sys = ss(A,B,C,0)
% rlocus(sys,linspace(-1,1))
K_P = .7;
eig(A+B*K_P*C)

% Choose a gain parameter eps, sufficiently small (test)
eps = .4;

Ae0 = [A+B*K_P*C,zeros(2,1);C,0];
Be0 = [B;0];
Ce0 = [zeros(1,2),1/(C*((-A)\B))];

sys = ss(Ae0,Be0,Ce0,0)
rlocus(sys,linspace(0,.27,3000))
%%

eps = 0.27;

[Ae,Be,Ce,De] = LinSysPIClosedLoop(A,B,C,K_P,eps);

% LinSysPlotEigs(Ae,[-1,0,-3,3])
% %


yref = @(t) -4;
yref = @(t) (-4)*(t<30) + (-2)*(t>=30);
% yref = @(t) 0.1*sin(0.1*t);


% The closed-loop system can be simulated with 'LinSysSim', now with the
% input function 'yref(t)'

% Initial state of the oscillator
x0 = [1;0];
% Initial state of the PI-controller
xc0 = 0;

tspan = [0 60];

sol = LinSysSim(Ae,Be,[x0;xc0],yref,tspan);

tt = linspace(tspan(1),tspan(2),500);
xxe = deval(sol,tt);

% The output of the controlled system is C*x(t) = [C,zeros(p)]*x_e(t)
yy = [C,0]*xxe;

% Values of yref(t) for plotting
yrefvals = zeros(1,length(tt));
for ind = 1:length(tt), yrefvals(ind)=yref(tt(ind)); end

figure(1)
% Plot the output and the reference
plot(tt,[yrefvals;yy],'Linewidth',2)
title(['Output for $K_P= ' num2str(K_P) '$ and $\varepsilon= ' num2str(eps) '$'],'Interpreter','Latex','Fontsize',16)

%% Animate the motion of the oscillator


figure(2)
massfig = [0,0,2,2;-1/6,1/6,1/6,-1/6];
axis([-6,6,-1,1])
hold on
for ind = 1:length(tt)
  cla
  
  patch(yy(ind)+massfig(1,:),massfig(2,:),[.7,0,0.5])
  plot(yy(ind)+[0,0],[0,-0.4],'k','linewidth',2)
  title(['time = ' num2str(round(tt(ind),1))])
  drawnow
  pause(0.01)
end

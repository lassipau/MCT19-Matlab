A = [0 1;-1 0];
B = [0;1];
C = [1 1/2];
D = 0;

x0 = [1;1];
tspan = [0 15];

ufun = @(t) sin(t).*cos(t);

sol = LinSysSim(A,B,x0,ufun,tspan);

%figure(1)
%LinSysStatePlot(sol,100,[tspan 1.1*[min(min(sol.y)) max(max(sol.y))]],2);
LinSysStatePlot(sol,100,[],2);
%figure(2)
%LinSysOutputPlot(sol,C,D,ufun,200);
%LinSysFigAdjust(axis)


%%


N = 8;
ga = -1;

A = ga*eye(N);
B = eye(N);
C = eye(N);
D = zeros(N);

x0 = exp(1i*((0:N-1)'*2*pi/N));

tspan = [0 6];

ufun = @(t) zeros(N,1);

sol = LinSysSim(A,B,x0,ufun,tspan);

%figure(1)
%LinSysStatePlot(sol,100,[tspan 1.1*[min(min(sol.y)) max(max(sol.y))]],2);
%LinSysStatePlot(sol,100,[],2);

tt = linspace(sol.x(1),sol.x(end),200);
xx = deval(sol,tt);

hold off
cla
hold on
plot(real(xx)',imag(xx)','color',.5*[1 1 1],'linewidth',2);
plot(real(xx(:,1))',imag(xx(:,1))','bo','markersize',10,'markerfacecolor',.5*[1 1 1],'markeredgecolor','k','linewidth',2);
set(gca,'box','on','tickdir','out','xtick',-1:1,'ytick',-1:1,'xticklabel',[],'yticklabel',[])
axis(1.1*[-1 1 -1 1])
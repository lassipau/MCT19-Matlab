N = 9;
ga = -1;

A = ga*eye(N);
B = eye(N);
C = eye(N);
D = zeros(N);

x0 = exp(1i*((0:N-1)'*2*pi/N));
%x0 = (1/2+1/2*rand(N,1)).*exp(1i*rand(N,1)*2*pi);

tspan = [0 6];

ufun = @(t) zeros(N,1);

sol = LinSysSim(A,B,x0,ufun,tspan);


tt = linspace(sol.x(1),sol.x(end),200);
xx = deval(sol,tt);


hold off
cla
hold on
plot(real(xx)',imag(xx)','color',.5*[1 1 1],'linewidth',2);
plot(real(xx(:,1))',imag(xx(:,1))','bo','markersize',10,'markerfacecolor',.5*[1 1 1],'markeredgecolor','k','linewidth',2);
axis off
axis(1.1*[-1 1 -1 1])
axis equal
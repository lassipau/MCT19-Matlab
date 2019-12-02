N = 9;
ga = 0;

A = ga*eye(N);
B = eye(N);
C = eye(N);
D = zeros(N);

x0 = exp(1i*((0:N-1)'*2*pi/N))
% x0 = (1/2+1/2*rand(N,1)).*exp(1i*rand(N,1)*2*pi);

% Stabilization and transfer function P_{K_P}(0)
% (system already stable, no need for additional stabilization
K_P = -eye(N);
PK0 = C*((-A-B*K_P*C)\B);

epsval = 0.3;
K_I = -epsval*(PK0\eye(N));

% Form the closed-loop system 
Ae = [A+B*K_P*C, B*K_I;C,zeros(N)];
Be = [-B*K_P;-eye(N)];
Ce = [C,zeros(N)];
De = -eye(N);

evals = eig(Ae);
if max(real(evals))>=0
  error('The closed-loop system is not stable! Adjust parameter epsval!')
end

% Initial state of the closed-loop system
% xe0 = [x0;zeros(N,1)];


tspan = [0 15];
randvec = rand(N,1).*exp(1i*rand(N,1)*2*pi);
yref = @(t) ones(N,1);
yref = @(t) x0;
yref = @(t) 2*exp(0.3*1i*t)*x0;
% yref = @(t) 2*exp(0.3*1i*t)*randvec;
% yref = @(t) 2*x0;
% yref = @(t) 2*randvec;
% yref = @(t) 2*exp(pi/2*1i)*x0;


xe0 = [x0;zeros(N,1)];
% xe0 = [zeros(N,1);zeros(N,1)];
% An adjusted initial state
xe0 = -Ae\Be*yref(0);


sol = LinSysSim(Ae,Be,xe0,yref,tspan);


tt = linspace(sol.x(1),sol.x(end),200);
xx = deval(sol,tt);

% plot(tt,xx((N+1):end,:))
xx = xx(1:N,:);

yrefplot = zeros(N,200);
% for ind = 1:200, yrefplot(:,ind) = yref(tt(ind)); end
% yrefplot = yref(tt);
% plot(real(yrefplot),imag(yrefplot))


% LinSysOutputPlot(sol,Ce,De,yref,200);

% sol.y = sol.y((N+1):end,:);
% LinSysStatePlot(sol,200)
% %
% hold off
% cla
% hold on
% plot(real(xx)',imag(xx)','color',.5*[1 1 1],'linewidth',2);
% plot(real(xx(:,1))',imag(xx(:,1))','bo','markersize',10,'markerfacecolor',.5*[1 1 1],'markeredgecolor','k','linewidth',2);
% axis off
% axis(3.1*[-1 1 -1 1])
% axis equal


%% Animate

tt = linspace(sol.x(1),sol.x(end),400);
xx = deval(sol,tt);

% plot(tt,xx((N+1):end,:))
xx = xx(1:N,:);

figure(2)
for ind = 1:length(tt)
  plot(real(xx(:,ind))',imag(xx(:,ind))','bo','markersize',10,'markerfacecolor',.5*[1 1 1],'markeredgecolor','k','linewidth',2);
  axis off
  axis(3.1*[-1 1 -1 1])
  axis equal
  title(['time = ' num2str(round(tt(ind),1))])
  drawnow
  pause(0.005)
end

%% plot the tracking error w.r.t time

tt = linspace(sol.x(1),sol.x(end),400);
xx = deval(sol,tt);
% y = C*xx;
error = zeros(N,length(tt));
errornorm = zeros(1,length(tt));
for ind = 1:length(tt)
  error(:,ind) = C*xx(1:N,ind)-yref(tt(ind)); 
  errornorm(ind) = norm(error(:,ind)); 
end

plot(tt,abs(error),'Linewidth',2)
% plot(tt,errornorm,'Linewidth',2)
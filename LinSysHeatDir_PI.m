% PI Control of a 1D Heat equation with Dirichlet boundary conditions. 
% Numerical approximation using "Finite Differences". 
% 
% The script constructs a PI controller for the heat equation
% 
% v_t(\xi,t) = \alpha*v_{\xi\xi}(\xi,t) + B(\xi)u(t)
% y(t) = \int_0^1 v(\xi,t)c(\xi)d\xi
%
% on [0,1] with Dirichlet boundary conditions v(0,t)=v(1,t)=0. 
% The approximation has the form x'(t)=A_N*x(t), where state of x(t) the 
% is a vector of the values (v(h,t),v(2*h,t),v(3*h,t), ..., v(1-h,t))^T for
% each time t. The length of the vector is N, and h=1/(N+1). The values at
% the boundaries v(0,t) and v(1,t) are not chosen as variables, since they
% are always zero.


% Size of the numerical approximation
N = 35;
% The heat conductivity of the material
alpha = .3;

% Define the system parameters. The approximation works best if all of the
% functions satisfy the boundary conditions, but you can still experiment
% functions that are not continuously differentiable, for example 
% f(x) = 1/2-abs(x-1/2)

% Define the input profile function b(x) and the output weight function
% c(x)
% b_fun = @(xi) 4*(xi>=1/4).*(xi<=1/2);
% c_fun = @(xi) 4*(xi>=1/2).*(xi<=3/4);

b_fun = @(xi) 10*(xi>=.15).*(xi<=.25);
c_fun = @(xi) 10*(xi>=.75).*(xi<=.85);


% Set the simulation parameters

% Initial state of the plant. You can define your initial conditions here!
% The input parameter 'xi' is the spatial variable (this should accept a 
% vector input)
x0fun = @(xi) zeros(size(xi));
x0fun = @(xi) 3*xi.*(1-xi);
% x0fun = @(xi) sin(pi*xi)+3*sin(5*pi*xi);


% Choose the simulation interval, 
tspan = [0,15];
% Choose the numbers of points in t- and x-directions for plotting
Nt = 200;


% Tracking: Reference output level, defined as a function of t
yref = @(t) -2*ones(size(t));
% yref = @(t) (-4)*(t<15) + (-2)*(t>=15);
% yref = @(t) 2*(t<15) + 4*(t>=15);
% yref = @(t) 0.3*sin(0.05*t)

% PI controller parameters

% The system is already exponentially stable, we can choose K_P=0
K_P = 0;

% If known exactly, the parameter P_{K_P}(0) can be entered here, 
% or it can be computed based on the numerical approximation later 
% (the latter is somewhat dangerous! But this does typically work for 
% heat equations.)
% PK0 = ...

% The value of epsilon, can use trial-and-error
epsgain = .5;





%%%%%%%%%%%%%%%%% Construction of the approximation %%%%%%%%%%%%%%%%%%%%


% The heat equation has Dirichlet boundary condition at x=0 and x=1. 
% The values of the solution at these two points are not chosen as variables

spgridfull = linspace(0,1,N+2);
h = spgridfull(2)-spgridfull(1);
spgrid = spgridfull(2:(end-1));

ee = ones(N,1);
AN = alpha*1/h^2*spdiags([ee -2*ee ee],-1:1,N,N);

% A perturbation:
gammafun = @(x) 5*(-abs(x-1/2)+1/2);
% AN = AN+spdiags(gammafun(spgrid).',0,N,N);
% plot(gammafun(0:0.1:1))
% %

% Convert the functions b(.) and c(.) to the matrices B and C of the 
% approximate system
BN = b_fun(spgrid.');
CN = h*c_fun(spgrid); % Integration over [0,1] is approximated with a sum


% Convert the initial condition to the initial state of the numerical
% approximation
x0 = x0fun(spgridfull(2:(N+1))).';

%%%%%%%%%%%%%%%%% Construct the PI-controller %%%%%%%%%%%%%%%%%%%%

% Compute PK0 based on the numerical approximation. Should always use
% DIFFERENT numerical approximations for computing P_{K_P}(0) and the
% simulation!
PK0 = -CN*((AN+BN*K_P*CN)\BN);
% Explicit value for control on [0.15,0.25] and observation on [0.75,0.85]
% PK0 = 0.04/alpha;

[Ae,Be,Ce,De] = LinSysPIClosedLoopInfDim(AN,BN,CN,K_P,PK0,epsgain);

% closed-loop initial state
xe0 = [x0;zeros(size(CN,1))];

%%%%%%%%%%%%%%%%% Simulation of the Closed-Loop System %%%%%%%%%%%%%%%%%%%%


% Complete the simulation. 'ode15s' solves a differential equation
% numerically
odefun = @(t,xe) Ae*xe + Be*yref(t);
sol = ode15s(odefun,tspan,xe0);
tgrid = linspace(tspan(1),tspan(2),Nt);
xxe = deval(sol,tgrid);


% Ae_evals = eig(full(Ae));
% set(gca,'xlim',[-4,0])
% 
% plot(real(Ae_evals),imag(Ae_evals),'b.','Markersize',20)
% %

%%%%%%%%%%%%%%%%% Illustration of the results %%%%%%%%%%%%%%%%%%%%


% Compute the output y(t) at the points 'tgrid' and plot the result
yy = Ce*xxe;

% Values of yref(t) for plotting
yrefvals = zeros(1,length(tgrid));
for ind = 1:length(tgrid), yrefvals(ind)=yref(tgrid(ind)); end

% Plot the output and the reference
figure(1)
plot(tgrid,[yrefvals;yy],'Linewidth',2)
title(['Output for $K_P= ' num2str(K_P) '$ and $\varepsilon= ' num2str(epsgain) '$'],'Interpreter','Latex','Fontsize',16)

%% Plot the solution of the heat equation as a function of t and xi

% To plot the solution on [0,1], we add the zero boundary values of the 
% solution at \xi=0 and \xi=1 to the approximate solution of the heat 
% equation (in the first N rows of the matrix 'xxe')
figure(2)
colormap jet
LinSysPlot1DSurf([zeros(1,length(tgrid));xxe(1:N,:);zeros(1,length(tgrid))],spgridfull,tgrid)


%%
% Animate the solution of the heat equation 
figure(3)
% No movie recording
[~,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xxe(1:N,:);zeros(1,length(tgrid))],spgridfull,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgridfull,tgrid,0.03,1);

% Numerical approximation of a 1D Heat equation with Dirichlet
% boundary conditions using "Finite Differences". Includes control input
% u(t) and measured output y(t).
% 
% The script defines a Finite Difference approximation of the heat equation
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
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

% Size of the numerical approximation
N = 25;
% The heat conductivity of the material
alpha = .3;

% Define the system parameters. The approximation works best if all of the
% functions satisfy the boundary conditions, but you can still experiment
% functions that are not continuously differentiable, for example 
% f(x) = 1/2-abs(x-1/2)

% Define the input profile function b(x) and the output weight function
% c(x)
b_fun = @(xi) 4*(xi>=1/4).*(xi<=1/2);
b_fun = @(xi) 4*(xi>=0).*(xi<=1/4);
% b_fun = @(xi) 1-xi;
c_fun = @(xi) 4*(xi>=1/2).*(xi<=3/4);


% Set the simulation parameters

% Initial state of the plant. You can define your initial conditions here!
% The input parameter 'xi' is the spatial variable (this should accept a 
% vector input)
x0fun = @(xi) zeros(size(xi));
x0fun = @(xi) 3*xi.*(1-xi);
% x0fun = @(xi) sin(pi*xi)+3*sin(5*pi*xi);


% Choose the simulation interval, 
tspan = [0,10];
% Choose the numbers of points in t- and x-directions for plotting
Nt = 80;

% Choose the input function
u_fun = @(t) zeros(size(t));
u_fun = @(t) ones(size(t));
u_fun = @(t) sin(t);




%%%%%%%%%%%%%%%%% Construction of the approximation %%%%%%%%%%%%%%%%%%%%


% The heat equation has Dirichlet boundary condition at x=0 and x=1. 
% The values of the solution at these two points are not chosen as variables

spgridfull = linspace(0,1,N+2);
h = spgridfull(2)-spgridfull(1);
spgrid = spgridfull(2:(end-1));

ee = ones(N,1);
AN = alpha*1/h^2*spdiags([ee -2*ee ee],-1:1,N,N);


% Convert the functions b(.) and c(.) to the matrices B and C of the 
% approximate system
B = b_fun(spgrid.');
C = h*c_fun(spgrid); % Integration over [0,1] is approximated with a sum


% Convert the initial condition to the initial state of the numerical
% approximation
x0 = x0fun(spgridfull(2:(N+1))).';



%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%


% Complete the simulation. 'ode15s' solves a differential equation
% numerically
odefun = @(t,x) AN*x + B*u_fun(t);
sol = ode15s(odefun,tspan,x0);
tgrid = linspace(tspan(1),tspan(2),Nt);
xx = deval(sol,tgrid);


%%%%%%%%%%%%%%%%% Illustration of the results %%%%%%%%%%%%%%%%%%%%


% Compute the output y(t) at the points 'tgrid' and plot the result
yy = C*xx;
figure(1)
plot(tgrid,yy,'Linewidth',2);
title('Output $y(t)$ of the heat equation','interpreter','latex','fontsize',16)

%% Plot the solution of the heat equation as a function of t and xi

% To plot the solution on [0,1], we add the zero boundary values of the 
% solution at \xi=0 and \xi=1 to the matrix 'xx' 
figure(2)
colormap jet
LinSysPlot1DSurf([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgridfull,tgrid)


%%
% Animate the solution of the heat equation 
figure(3)
% No movie recording
[~,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgridfull,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgridfull,tgrid,0.03,1);

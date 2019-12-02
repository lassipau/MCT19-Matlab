% Numerical approximation of a 1D Heat equation with Neumann
% boundary conditions using "Finite Differences".
% 
% The script defines a Finite Difference approximation of the heat equation
% 
% v_t(\xi,t) = \alpha*v_{\xi\xi}(\xi,t) + f(\xi,t)
%
% on [0,1] with Dirichlet boundary conditions v_\xi(0,t)=v_\xi(1,t)=0. 
% The approximation has the form x'(t)=A_N*x(t), where state of x(t) the 
% is a vector of the values (v(0,t),v(h,t),v(2*h,t), ..., v(1-h,t), v(1,t))^T for
% each time t. The length of the vector is N, and h=1/(N-1). 


% Size of the numerical approximation
N = 30;
% The heat conductivity of the material
alpha = .1;


% The heat equation has Dirichlet boundary condition at x=0 and x=1. 

spgrid = linspace(0,1,N);
h = spgrid(2)-spgrid(1);

ee = ones(N,1);
AN = alpha*1/h^2*spdiags([ee -2*ee ee],-1:1,N,N);
% Adjust the matrix to model the Neumann boundary conditions
AN(1,2) = alpha*2/h^2;
AN(end,end-1) = alpha*2/h^2;


% Initial state of the plant. You can define your initial conditions here!
% The input parameter 'xi' is the spatial variable (this should accept a 
% vector input)
% x0fun = @(xi) zeros(size(xi));
x0fun = @(xi) xi.*(1-xi);

% The "source function" of the heat equation, you can define your own! The
% input parameter 'xi' is the spatial variable (this should accept a vector 
% input) and t is the time.
% source_fun = @(xi,t) zeros(size(xi));
source_fun = @(xi,t) ones(size(xi))*sin(t);
% source_fun = @(xi,t) xi.*(1-xi).^2;
% source_fun = @(xi,t) (1/2-abs(xi-1/2))*t;
% source_fun = @(xi,t) 5*xi.*(1-xi).^3.*(t<5);
% source_fun = @(xi,t) 5*xi.*(1-xi).^3.*(rem(t,2)<1);


% Convert the initial condition to the initial state of the numerical
% approximation
x0 = x0fun(spgrid).';

% Simulate the approximation of the heat equation

% Define the time-interval of the simulation
tspan = [0,20];


odefun = @(t,x) AN*x + source_fun(spgrid.',t);
sol = ode15s(odefun,tspan,x0);
tgrid = linspace(tspan(1),tspan(2),100);
xx = deval(sol,tgrid);

% Plot the solution of the heat equation as a function of t and xi
% To plot the solution on [0,1], we add the zero boundary values of the 
% solution at \xi=0 and \xi=1 to the matrix 'xx' 
figure(1)
colormap jet
LinSysPlot1DHeatSurf(xx,spgrid,tgrid)


%%
% Animate the solution of the heat equation 
figure(2)
% To animate the solution on [0,1], we add the zero boundary values of the 
% solution at \xi=0 and \xi=1 to the matrix 'xx' 
% No movie recording
[~,zlims] = LinSysAnimate1D(xx,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgrid,tgrid,0.03,1);

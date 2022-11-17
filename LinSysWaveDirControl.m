% Numerical approximation of a 1D Wave equation with Dirichlet
% boundary conditions using "Modal Approximation". Includes control input
% u(t) (force input) and measured output y(t) (weighted average of
% velocity).
% 
% The script defines a Modal approximation of the heat equation
% 
% w_{tt}(\xi,t) = c^2*w_{\xi\xi}(\xi,t) - d(\xi)*w_t(\xi,t) + B(\xi)u(t)
% y(t) = \int_0^1 w_\xi(\xi,t)*c_1(\xi)d\xi+\int_0^1 w_t(\xi,t)*c_2(\xi)d\xi
%
% on [0,1] with Dirichlet boundary conditions v(0,t)=v(1,t)=0. 
% The approximation has the form x'(t)=A_N*x(t)+B_N*u(t), where state of x(t) the 
% is a vector of the coordinates (\alpha_1(t),...,\alpha_N(t))^T of the
% approximate wave profile in the basis spanned by the eigenvectors of the
% second order differential operator. For plotting the solution the values
% of these linear combinations are evaluated at a spatial grid on [0,1].
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)


% Size of the numerical approximation
N = 45;
% The wave speed
c = .3;

% Define the system parameters. The approximation works best if all of the
% functions satisfy the boundary conditions, but you can still experiment
% functions that are not continuously differentiable, for example 
% f(xi) = 1/2-abs(xi-1/2)
% Define the damping profile function
d = @(xi) 5*xi.*(1-xi);
% d = @(xi) 1/2-abs(xi-1/2);
% d = @(xi) zeros(size(xi));
% d = @(xi) 0.5*ones(size(xi));

% Define the input profile function b(xi) and the output weight function c_2
b = @(xi) 10*xi.^2.*(1-xi);
c_1 = @(xi) xi.*(1-xi);
c_1 = @(xi) 50*(xi<=0.11).*(xi>=0.09);

c_2 = @(xi) xi.*(1-xi);
% c_2 = @(xi) zeros(size(xi));


% Set the simulation parameters

% Choose the initial wave profile w_0(\xi) and initial velocity w_1(\xi)
w_0 = @(xi) sin(pi*xi)+3*sin(5*pi*xi);
w_0 = @(xi) 10*xi.*(1-xi)-2*sin(2*pi*xi);
% w_0 = @(xi) 1/2-abs(xi-1/2);

w_1 = @(xi) zeros(size(xi));

% Choose the simulation interval, 
tspan = [0,15];
% Choose the numbers of points in t- and x-directions for plotting
Nt = 200;
Nx = 160;

% Choose the input function
u_fun = @(t) zeros(size(t));
% u_fun = @(t) ones(size(t));
% u_fun = @(t) sin(t);


%%%%%%%%%%%%%%%%% Construction of the approximation %%%%%%%%%%%%%%%%%%%%

% Eigenvalues of the second order differential operator (d^2/d\xi^2) with
% Dirichlet boundary conditions
lambda_n = @(n) -c^2*n.^2*pi^2;
% Orthonormal eigenvectors of the operator (d^2/d\xi^2) with Dirichlet 
% boundary conditions
phi_n = @(xi,n) sqrt(2)*sin(n*pi*xi);
phi_n_prime = @(xi,n) -sqrt(2)*n*pi*cos(n*pi*xi);


% Define a Matlab function for computing the inner product <f,\phi_n>,
% which gives the n'th coordinate of the function f(x) in the orthonormal 
% basis spanned by the eigenvectors \phi_n(x) (this is the Fourier sine
% series).
phiprod = @(f,n) quadgk(@(xi) f(xi).*phi_n(xi,n),0,1,'AbsTol',1e-9);
phiprod_prime = @(f,n) quadgk(@(xi) f(xi).*phi_n_prime(xi,n),0,1,'AbsTol',1e-9);


% Construct the matrices A_N, B_N and C_N of the approximate system
A_0 = spdiags(lambda_n((1:N).'),0,N,N);

B_0 = zeros(N,1);
C_1 = zeros(1,N);
C_2 = zeros(1,N);
for ind = 1:N
  B_0(ind) = phiprod(b,ind);
  C_1(ind) = phiprod_prime(c_1,ind);
  C_2(ind) = phiprod(c_2,ind);
end

D_0 = zeros(N,N);
for ind1 = 1:N
  for ind2 = ind1:N
    D_0(ind1,ind2) = phiprod(@(xi) d(xi).*phi_n(xi,ind2),ind1);
  end
  D_0((ind1+1):end,ind1) = D_0(ind1,(ind1+1):end).';
end


AN = [sparse(N,N),speye(N);A_0,-D_0];
BN = [zeros(N,1);B_0];
CN = [C_1,C_2];


%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%


% Define the initial state of the approximate system. 
w0_vec = zeros(N,1);
w1_vec = zeros(N,1);
for ind = 1:N
  w0_vec(ind) = phiprod(w_0,ind);
  w1_vec(ind) = phiprod(w_1,ind);
end
alpha0 = [w0_vec;w1_vec];


% Complete the simulation. 'ode15s' solves a differential equation
% numerically
odefun = @(t,alpha) AN*alpha + BN*u_fun(t);
sol = ode15s(odefun,tspan,alpha0);
tgrid = linspace(tspan(1),tspan(2),Nt);
alphas = deval(sol,tgrid);



%%%%%%%%%%%%%%%%% Illustration of the results %%%%%%%%%%%%%%%%%%%%

% Compute the output y(t) at the points 'tgrid' and plot the result
yy = CN*alphas;
figure(1)
plot(tgrid,yy,'Linewidth',2);
title('Output $y(t)$ of the heat equation','interpreter','latex','fontsize',16)

%% Plot the solution of the heat equation as a function of t and xi

% Compute the values of the approximate solution from the coordinates in
% the matrix 'alphas' (the coordinates for the deflection w(\xi,t) are in 
% the N first rows). This can be achieved by multiplying N first rows of 
% 'alphas' from the left with a suitable matrix 'PHI' whose columns contain
% the values of each of the basis functions 'phi_n(\xi)' in the spatial 
% grid on [0,1].
% Compute the matrix 'PHI'
spgrid = linspace(0,1,Nx);
PHI = zeros(Nx,N);
for ind = 1:N
  PHI(:,ind) = phi_n(spgrid.',ind);
end
zz = PHI*alphas(1:N,:);

% To plot the solution on [0,1], we add the zero boundary values of the 
% solution at \xi=0 and \xi=1 to the matrix 'xx' 
figure(2)
colormap jet
LinSysPlot1DSurf(zz,spgrid,tgrid)


%%
% Animate the solution of the heat equation 
figure(3)
% No movie recording
[~,zlims] = LinSysAnimate1D(zz,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = LinSysAnimate1D([zeros(1,length(tgrid));xx;zeros(1,length(tgrid))],spgridfull,tgrid,0.03,1);

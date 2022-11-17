% Numerical approximation of a 1D Wave equation with Dirichlet
% boundary conditions using "Modal Approximation". Includes control input
% u(t) (force input) and measured output y(t) (weighted average of
% velocity).
% 
% The script defines a Modal approximation of the heat equation
% 
% w_{tt}(\xi,t) = c^2*w_{\xi\xi}(\xi,t) - d(\xi)*w_t(\xi,t) + B(\xi)u(t)
% y(t) = \int_0^1 w(\xi,t)*c_1(\xi)d\xi
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
% NOTE THAT if you change these, you also need to change the value of PK0
% computed later in the code (around line 150).
% b = @(xi) 10*xi.^2.*(1-xi);
b = @(xi) 10*(xi>=0.1).*(xi<=0.3).*(1-10*abs(xi-0.2));

% c_1 = @(xi) xi.*(1-xi);
% c_1 = @(xi) 50*(xi<=0.11).*(xi>=0.09);
c_1 = @(xi) 10*(xi>=0.7).*(xi<=0.9).*(1-10*abs(xi-0.8));

% Plot the input and output functions
% xx = linspace(0,1,401);
% plot(xx,[b(xx);c_1(xx)],'Linewidth',2);
% axis([0,1,0,11])
% %

% Set the simulation parameters

% Choose the initial wave profile w_0(\xi) and initial velocity w_1(\xi)
% w_0 = @(xi) sin(pi*xi)+3*sin(5*pi*xi);
w_0 = @(xi) 10*xi.*(1-xi)-2*sin(2*pi*xi);
% w_0 = @(xi) 1/2-abs(xi-1/2);
% w_0 = @(xi) zeros(size(xi));

w_1 = @(xi) zeros(size(xi));

% Choose the simulation interval, 
tspan = [0,30];
% Choose the numbers of points in t- and x-directions for plotting
Nt = 200;
Nx = 160;

% Better settings for surface plots
% Nt = 100;
% Nx = 60;

% Tracking: Reference output level, defined as a function of t
yref = @(t) 1*ones(size(t));
% yref = @(t) (-4)*(t<15) + (-2)*(t>=15);
% yref = @(t) 2*(t<15) + 4*(t>=15);
% yref = @(t) 0.3*sin(0.05*t)

% PI controller parameters

% The system is already exponentially stable, we can choose K_P=0
K_P = 0;

% If known exactly, the parameter P_{K_P}(0) can be entered in the part
% "Construct the PI-controller" below. Alternatively it can be computed 
% using a numerical approximation. The latter is quite dangerous in the 
% case of wave equations, since these PDEs do not have good approximation 
% properties! 

% The value of epsilon, can use trial-and-error ("best" around .25)
epsgain = .15;
epsgain = .25;
% epsgain = .35;


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
for ind = 1:N
  B_0(ind) = phiprod(b,ind);
  C_1(ind) = phiprod(c_1,ind);
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
CN = [C_1,zeros(1,N)];


% Define the initial state of the approximate system. 
w0_vec = zeros(N,1);
w1_vec = zeros(N,1);
for ind = 1:N
  w0_vec(ind) = phiprod(w_0,ind);
  w1_vec(ind) = phiprod(w_1,ind);
end
alpha0 = [w0_vec;w1_vec];



%%%%%%%%%%%%%%%%% Construct the PI-controller %%%%%%%%%%%%%%%%%%%%

% Compute PK0 based on the numerical approximation. If this is done 
% numerically, you should always use a DIFFERENT numerical approximations 
% for computing P_{K_P}(0) and the simulation!
% PK0 = -CN*((AN+BN*K_P*CN)\BN);
% Explicit value for control in Example 5.2.3
PK0 = 1/(25*c^2);

[Ae,Be,Ce,De] = LinSysPIClosedLoopInfDim(AN,BN,CN,K_P,PK0,epsgain);

% closed-loop initial state
xe0 = [alpha0;zeros(size(CN,1))];

% Ae_evals = eig(full(Ae));
% set(gca,'xlim',[-4,0])
% 
% plot(real(Ae_evals),imag(Ae_evals),'b.','Markersize',20)
% %

%%%%%%%%%%%%%%%%% Simulation of the Closed-Loop System %%%%%%%%%%%%%%%%%%%%




% Complete the simulation. 'ode15s' solves a differential equation
% numerically
odefun = @(t,xe) Ae*xe + Be*yref(t);
sol = ode15s(odefun,tspan,xe0);
tgrid = linspace(tspan(1),tspan(2),Nt);
xxe = deval(sol,tgrid);

% Extract the state of the controlled system, in columns 1..2*N of 'xxe'
alphas = xxe(1:(2*N),:);


%%%%%%%%%%%%%%%%% Illustration of the results %%%%%%%%%%%%%%%%%%%%

% Compute the output y(t) at the points 'tgrid' and plot the result
yy = CN*alphas;

% Values of yref(t) for plotting
yrefvals = zeros(1,length(tgrid));
for ind = 1:length(tgrid), yrefvals(ind)=yref(tgrid(ind)); end

% Plot the output and the reference
figure(1)
plot(tgrid,[yrefvals;yy],'Linewidth',2)
title(['Output for $K_P= ' num2str(K_P) '$ and $\varepsilon= ' num2str(epsgain) '$'],'Interpreter','Latex','Fontsize',16)

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

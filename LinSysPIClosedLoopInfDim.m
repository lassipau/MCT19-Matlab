function [Ae,Be,Ce,De] = LinSysPIClosedLoopInfDim(A,B,C,K_P,PK0,epsgain)
% function [A_e,B_e,C_e,D_e] = LinSysPIClosedLoop(A,B,C,K_P,PK0,epsgain)
% 
% Form the closed-loop system (Ae,Be,Ce,De) consisting of the linear system
% (A,B,C) and a Proportional-Integral Controller (PI Controller) with the
% parameters K_P (proportional part gain) K_I = eps*(C*(A+B*K_P*C)^{-1}B)^{-1} 
% (integral part gain) where eps>0 is a low-gain parameter. The routine
% tests the stability of the closed-loop system.
%
% Parameters: 
% A = nxn-matrix, B = nxm-matrix, C = pxn-matrix, 
% K_P = mxp-matrix, eps>0
% PK0 = pxm-matrix (an approximate) value of P_{K_P}(0)
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

p = size(C,1);
m = size(B,2);

if ~isequal(size(K_P),[m,p])
    error('K_P has incorrect dimensions!')
elseif ~isequal(size(PK0),[p,m])
    error('PK0 has incorrect dimensions!')
end

if find(real(eig(A+B*K_P*C))>=0)
  warning('The matrix A+B*K_P*C is not Hurwitz!')
end

if rank(PK0,1e-10)<p
  error('The transfer function of (A,B,C) is nearly non-surjective at s=0!')
end

K_I = -epsgain*pinv(PK0);

Ae = [A+B*K_P*C,B*K_I;C,zeros(p)];
Be = [-B*K_P;-eye(p)];
Ce = [C,zeros(p)];
De = -eye(p);

% Test the stability of the closed-loop system, and print out the stability
% margin.
CLeigs = eig(Ae);
maxRe = max(real(CLeigs));

if maxRe>=0
  error('The closed-loop system matrix Ae is not stable! Adjust controller parameters!')
end

fprintf(['The largest real part of eigenvalues of Ae = ' num2str(maxRe) '\n'])


function [A_K,B_K,C_K,D_K] = LinSysOutputFB(A,B,C,D,K)
% function [A_K,B_K,C_K,D_K] = LinSysOutputFB(A,B,C,D,K)
%
% Copyright (C) 2019 by Lassi Paunonen (lassi.paunonen@tuni.fi)

p = size(C,1);
IDK = eye(p)-D*K;

if rank(IDK)<n
    error('The matrix I-DK is not invertible')
end

A_K = A + B*K*(IDK\C);
B_K = B/(eye(size(B,2)-K*D));
C_K = IDK\C;
D_K = IDK\D;
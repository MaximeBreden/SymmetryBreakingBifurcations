function [F, DuF, DlambdaF, DlambdauFphi, DuuFphi] = F_DF_SKT(U, para, N, Phi)

% The map F corresponding to the SKT system, and some of its derivatives.

%% Initialization

% Check whether we are using Intlab and intervals or not
if isa(U(1), 'intval')
    ipi = intval('pi');
else
    ipi = pi;
end

% check whether the solution is 1D or 2D
is1D = not(isfield(para,'Ly')); 
    
NN = prod(N);
U1 = U(1:NN);
U2 = U(NN+1:end);

d1 = para.d1;
d2 = para.d2;
d11 = para.d11;
d12 = para.d12;
d21 = para.d21;
d22 = para.d22;
r1 = para.r1;
a1 = para.a1;
b1 = para.b1;
r2 = para.r2;
a2 = para.a2;
b2 = para.b2;
L1 = para.Lx;
if not(is1D)
    L2 = para.Ly;
end

%% F
Mu1 = convomat(reshape(U1, N));
Mu2 = convomat(reshape(U2, N));
U1U1 = Mu1*U1;
U2U2 = Mu2*U2;
U1U2 = Mu1*U2;
if is1D
    Lap = -(ipi/L1)^2 * ((0:N(1)-1)').^2;
else
    [N2, N1] = meshgrid(0:N(2)-1,0:N(1)-1);
    Lap = -(ipi/L1)^2*N1.^2 - (ipi/L2)^2*N2.^2;
    Lap = Lap(:);
end

Dif1 = d1*U1 + d11*U1U1 + d12*U1U2;
Dif2 = d2*U2 + d21*U1U2 + d22*U2U2;
R1 = r1*U1 - a1*U1U1 - b1*U1U2;
R2 = r2*U2 - b2*U1U2 - a2*U2U2;

F = [Lap.*Dif1 + R1;
     Lap.*Dif2 + R2];

%% DuF
if nargout > 1 %if we also want DuF 
    I = eye(NN);  

    DDif11 = d1*I + 2*d11*Mu1 + d12*Mu2;
    DDif12 = d12*Mu1;
    DDif21 = d21*Mu2;
    DDif22 = d2*I + d21*Mu1 + 2*d22*Mu2;

    DR11 = r1*I - 2*a1*Mu1 - b1*Mu2;
    DR12 = -b1*Mu1;
    DR21 = -b2*Mu2;
    DR22 = r2*I - b2*Mu1 - 2*a2*Mu2;

    DuF = [Lap.*DDif11 + DR11, Lap.*DDif12 + DR12;
           Lap.*DDif21 + DR21, Lap.*DDif22 + DR22];
end

%% Extra derivatives
if nargout > 2 %if we need extra derivatives for the extended system
    if nargin < 3
        error("You need to input the eigenfunction phi if you want these derivatives")
    end
    Phi1 = Phi(1:NN);
    Phi2 = Phi(NN+1:end);

    DlambdaF = [Lap.*U1;
                Lap.*U2];

    DlambdauFphi = [Lap.*Phi1;
                    Lap.*Phi2];

    Mphi1 = convomat(reshape(Phi1, N));
    Mphi2 = convomat(reshape(Phi2, N));
    DuuFphi = [Lap.*(2*d11*Mphi1+d12*Mphi2)-2*a1*Mphi1-b1*Mphi2, Lap.*(d12*Mphi1)-b1*Mphi1;
               Lap.*(d21*Mphi2)-b2*Mphi2, Lap.*(d21*Mphi1+2*d22*Mphi2)-b2*Mphi1-2*a2*Mphi2];
end


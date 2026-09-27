function [F, DuF, DlambdaF, DlambdauFphi, DuuFphi] = F_DF_OK(U, para, N, Phi)

% The map F corresponding to the Ohta-Kawasaki equation, and some of its 
% derivatives.

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
u = reshape(U, N);

lambda = para.lambda;
mu = para.mu;
sigma = para.sigma;
L1 = para.Lx;
if not(is1D)
    L2 = para.Ly;
end

%% F
Mu = convomat(u);
UU = Mu*U;
UUU = Mu*UU;
if is1D
    Lap = -(ipi/L1)^2 * ((0:N(1)-1)').^2;
else
    [N2, N1] = meshgrid(0:N(2)-1,0:N(1)-1);
    Lap = -(ipi/L1)^2*N1.^2 - (ipi/L2)^2*N2.^2;
    Lap = Lap(:);
end
E = zeros(NN,1);
E(1) = 1;
F = (-Lap.^2).*U - lambda*(Lap.*(U-UUU)) - lambda*sigma*(U-mu*E);

%% DuF
if nargout > 1 
    I = eye(NN);  
    Muu = convomat(reshape(UU,N));
    DuF = -diag(Lap.^2) - lambda*Lap.*(I-3*Muu) - lambda*sigma*I;
end

if nargout > 2 %if we need extra derivatives for the extended system
    if nargin < 3
        error("You need to input the eigenfunction phi if you want these derivatives")
    end
    DlambdaF = -(Lap.*(U-UUU)) - sigma*(U-mu*E);
    DlambdauFphi = -Lap.*((I-3*Muu)*Phi) - sigma*Phi;
    DuuFphi = 6*lambda*Lap.*convomat(reshape(Mu*Phi,N));
end


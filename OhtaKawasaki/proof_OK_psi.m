function [rmin, Abar] = proof_OK_psi(Z, X, para, N, nsym, nu, s, rX, Abar)

% Tries to prove the eigenpair of D_uF(lambda,u)^* described by Z=(mu; psi),
% i.e., to prove that G has a zero near Z.
% The last argument is optional, and is computed inside this function if 
% not given as input.

%% Initialization
% Check whether we are using Intlab and intervals or not
if isa(X(1), 'intval')
    useintervals = true;
    ipi = intval('pi');
else
    useintervals = false;
    ipi = pi;
end

% check whether the solution is 1D or 2D
is1D = not(isfield(para,'Ly'));

% get lambda, u and phi
[lambda, u, phi] = extract_OK(X, N, nsym, 'mat');
U = u(:);
Phi = phi(:);

% This should already be true, but better be safe
para.lambda = lambda;

% get mu and psi
[mu, Psi] = extract_OK_psi(Z, N, 'vect');

% All the parameters of the system
sigma = para.sigma;
L1 = para.Lx;
if not(is1D)
    L2 = para.Ly;
end

% The approximate inverse A is defined by Abar
if nargin < 9
    [~, DuF] = F_DF_OK(U, para, N);
    Abar = inv([0, -Psi'; Phi, DuF-mu*eye(prod(N))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% We first derive all the bounds pretending that X is the exact solution,
%%% and then at the end we add all the errors terms coming from the fact 
%%% that X is only rX close to the exact solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Y
N3 = 3*N-2;
X3N = projection_OK(X, N3, N, nsym);
Z3N = projection_OK_psi(Z, N3, N);
G3N = G_DG_OK(Z3N, X3N, para, N3, nsym);
AadjG = applyAadj_OK_psi(Abar, G3N, para, N, N3);
Y = norm_blocks_OK_psi(AadjG, nu, s-4, N3);

% % Displaying the obtained estimates
% fprintf("\nY =\n")
% disp(Y)

%% Z1
%%% Finite part
[~, U3N, Phi3N] = extract_OK(X3N, N3, nsym, 'vectnosym');
[~, DuF3N] = F_DF_OK(U3N, para, N3);
ind_NinN3 = getindwithin(N, N3);
DuF3N = DuF3N(:,ind_NinN3); 

NN = prod(N);
NN3 = prod(N3);
Nfull = 1 + NN;
N3full = 1 + NN3;
I_3NN = zeros(N3full, Nfull);
I_3NN(1) = 1;
I_3NN(1+ind_NinN3,1+(1:NN)) = eye(NN);

L = [0, -Psi';
     Phi3N, DuF3N-mu*I_3NN(2:end,2:end)];

B = I_3NN - L*Abar;
Z1_finite = opnorm_blocks_OK_psi(B, nu, s-4, s-4, N3, N);
Z1_finite = Z1_finite';

% % Displaying the obtained estimates
% fprintf("\nZ1_finite =\n")
% disp(Z1_finite)

%%% tail part 
N2 = 2*N-1;
eN2 = zeros(N2);
eN2(1,1) = 1;

if s < 2
    error("The upcoming estimate assumes that \ell^1_{s-2} is a Banach algebra, and hence that s>=2.")
end
alpha = eN2(:) - 3*convomat(u, N2, N)*u(:);
norm_alpha = abs(alpha)' * getweights(N2, nu, s-2); % We use the s-2 norm here!
if is1D
    lambdaN = (ipi/L1)^2*N(1)^2;
    coef_alpha = ( N(1)/(N(1)-1) )^2;
    cstemin = (ipi/L1)^2;
else
    lambdaN = min((ipi/L1)^2*N(1)^2, (ipi/L2)^2*N(2)^2);
    mN =  min(N);
    coef_alpha = (max(L1,L2)/min(L1,L2))^2 * ((1+mN)/mN)^2;
    cstemin = (ipi/(min(L1,L2)))^2;
end

Z1_tail = zeros(2,2);
if useintervals
    Z1_tail = intval(Z1_tail);
end
Z1_tail(2,2) = abs(lambda)*norm_alpha*coef_alpha/lambdaN + abs(lambda*sigma+mu)/lambdaN^2;
Z1_tail = Z1_tail';

% % Displaying the obtained estimates
% fprintf("\nZ1_tail =\n")
% disp(Z1_tail)

%%% full Z1
Z1 = max(Z1_finite, Z1_tail);

% % Displaying the obtained estimates
% fprintf("\nZ1 =\n")
% disp(Z1)
% fprintf("\nSpectral radius of Z_1 : %f \n",eigs(i2f(Z1), 1)) % If this is >= 1, the proof will fail (it may also fail if this is too close to 1)
   
%% Z2 
D2G = zeros(2,2,2);
% the mu*psi term 
D2G(2,1,2) = 1;
D2G(2,2,1) = 1;

% norms of A
normA = opnorm_blocks_OK_psi(Abar, nu, s-4, s-4, N, N);
normA(2,2) = max(normA(2,2), 1/lambdaN^2);

Z2 = reshape( normA' * reshape(D2G, [2, 4]), [2, 2, 2]);

% % Displaying the obtained estimates
% fprintf("\nZ2 =\n")
% disp(Z2)

%% Error terms coming from the fact that G involves DuF at the exact solution

% Remember that the error bounds in rX are in the strong norm (\ell^1_s)
rlambda = rX(1);
ru = rX(2);

norm_u = abs(u(:))' * getweights(N, nu, s-2); % We use the s-2 norm here!

% The error terms (coming from DuF(X)-DuF(\bX))
errDFnoLap = [0, 0;
            0, rlambda*abs(sigma)];
errDFLap = [0, 0;
            0, cstemin * ( rlambda*norm_alpha + 3*(abs(lambda)+rlambda)*ru*(2*norm_u+ru) )];

% norm of A from k-4 to k-2
normLapA = opnorm_blocks_OK_psi(Abar, nu, s-2, s-4, N, N);
if is1D
    tailterm = (L1/ipi)^2 * ((1+N(1))/N(1))^2; 
else
    mN = min(N); 
    tailterm = (max(L1,L2)/ipi)^2 * ((1+mN)/mN)^2; 
end
normLapA(2,2) = max(normLapA(2,2), tailterm/lambdaN);


% DFA
errAadjepsadj = ( errDFnoLap*normA + errDFLap*normLapA )';

% Updated bounds
norm_psi = max(abs(Psi)./getweights(N, nu, s-4));
Yerr = errAadjepsadj * [0;norm_psi];
Y = Y + Yerr;

Z1err = errAadjepsadj;
Z1 = Z1 + Z1err;

%% Checking that we have a contraction
rmin = polynomialsnegative(Y,Z1,Z2);


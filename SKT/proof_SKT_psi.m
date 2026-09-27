function [rmin, Abar] = proof_SKT_psi(Z, X, para, N, nsym, nu, s, rX, w, Abar)

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
[lambda, u1, u2, phi1, phi2] = extract_SKT(X, N, nsym, 'mat');
U = [u1(:); u2(:)];
Phi = [phi1(:); phi2(:)];

% This should already be true, but better be safe
para.lambda = lambda;

% get mu and psi
[mu, Psi1, Psi2] = extract_SKT_psi(Z, N, 'vect');
Psi = [Psi1; Psi2];

% All the parameters of the system
d1 = lambda;
d2 = lambda;
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

% The approximate inverse A is defined by Abar
if nargin < 10
    [~, DuF] = F_DF_SKT(U, para, N);
    Abar = inv([0, -Psi'; Phi, DuF-mu*eye(2*prod(N))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% We first derive all the bounds pretending that X is the exact solution,
%%% and then at the end we add all the errors terms coming from the fact 
%%% that X is only rX close to the exact solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Y
N2 = 2*N-1;
NN2 = prod(N2);
N3 = 3*N-2;
NN3 = prod(N3);

X2N = projection_SKT(X, N2, N, nsym);
Z2N = projection_SKT_psi(Z, N2, N);
G2N = G_DG_SKT(Z2N, X2N, para, N2, nsym);
AadjG = applyAadj_SKT_psi(Abar, w, G2N, para, N, N2, N3);
Y = norm_blocks_SKT_psi(AadjG, nu, s-2, N3);

% % Displaying the obtained estimates
% fprintf("\nY =\n")
% disp(Y)

%% Z1
%%% Finite part
X3N = projection_SKT(X, N3, N, nsym);
[~, U13N, U23N, Phi13N, Phi23N] = extract_SKT(X3N, N3, nsym, 'vectnosym');
U3N = [U13N;U23N];
Phi3N = [Phi13N;Phi23N];

[~, Psi12N, Psi22N] = extract_SKT_psi(Z2N, N2, 'vect');
Psi2N = [Psi12N; Psi22N];

[~, DuF3N] = F_DF_SKT(U3N, para, N3);
ind_N2inN3 = getindwithin(N2, N3);
ind12_N2inN3 = [ind_N2inN3; NN3+ind_N2inN3];
DuF3N = DuF3N(:,ind12_N2inN3); 

N2full = 1 + 2*NN2;
N3full = 1 + 2*NN3;
I_3N2N = zeros(N3full, N2full);
I_3N2N(1) = 1;
I_3N2N(1+ind_N2inN3,1+(1:NN2)) = eye(NN2);
I_3N2N(1+NN3+ind_N2inN3,1+NN2+(1:NN2)) = eye(NN2);

Ladj = [0, Phi3N';
        -Psi2N, DuF3N'-mu*I_3N2N(2:end,2:end)'];

ind_NinN2 = getindwithin(N, N2);
inf_NinN2full = [1; 1+ind_NinN2; 1+NN2+ind_NinN2];
B = I_3N2N(:,inf_NinN2full) - applyAadj_SKT_psi(Abar, w, Ladj, para, N, N2, N)';
Z1_finite = opnorm_blocks_SKT_psi(B, nu, s-2, s-2, N3, N);
Z1_finite = Z1_finite';

% % Displaying the obtained estimates
% fprintf("\nZ1_finite =\n")
% disp(Z1_finite)

%%% tail part 
if s < 2
    error("The upcoming estimate assumes that \ell^1_{s-2} is a Banach algebra, and hence that s>=2.")
end
e = zeros(N);
e(1,1) = 1;
alpha = {d1*e+2*d11*u1+d12*u2, d12*u1;
         d21*u2, d2*e+d21*u1+2*d22*u2}; 
R = {r1*e-2*a1*u1-b1*u2, -b1*u1;
     -b2*u2, r2*e-b2*u1-2*a2*u2}; 

eN2 = zeros(N2);
eN2(1,1) = 1;
delta = myprod2x2(alpha,w,N2);
delta = {eN2-delta{1,1},-delta{1,2};
         -delta{2,1},eN2-delta{2,2}};

eta = myprod2x2(R, w, N2);

normdelta = norm2x2(delta, nu, s); % We use the s norm here!
normeta = norm2x2(eta, nu, s-2);
normw = norm2x2(w, nu, s-2);


Z1_tail = zeros(3,3);
if useintervals
    Z1_tail = intval(Z1_tail);
end
if is1D
    lambdaN = (ipi/L1)^2*N(1)^2;
    cstemax = (L1/ipi)^2;
    cstemin = (ipi/L1)^2;
    mN = N(1); % In the 1D case N is actually [N 1]
else
    lambdaN = min((ipi/L1)^2*N(1)^2, (ipi/L2)^2*N(2)^2);
    cstemax = (max(L1,L2)/ipi)^2;
    cstemin = (ipi/(min(L1,L2)))^2;
    mN = min(N);
end
Z1_tail(2:3,2:3) = cstemin*normdelta*cstemax*((1+mN)/mN)^2 + (normeta + abs(mu)*normw) / lambdaN;

psi1 = reshape(Psi1, N);
psi2 = reshape(Psi2, N);
q1 = convo(w{1,1}, psi1, N2) + convo(w{2,1}, psi2, N2);
q2 = convo(w{1,2}, psi1, N2) + convo(w{2,2}, psi2, N2);
if is1D
    Lap2N = (ipi/L1)^2 * ((0:N2(1)-1)').^2;
else
    [n2, n1] = meshgrid(0:N2(2)-1, 0:N2(1)-1);
    Lap2N = (ipi/L1)^2*n1.^2 + (ipi/L2)^2*n2.^2;
    Lap2N = Lap2N(:);
end
[~, ind_tail] = getindwithin(N, N2);
weights_Y2N = getweights(N2, nu, s-2);
Z1_tail(1,2) = max( abs(q1(ind_tail)) ./ (Lap2N(ind_tail) .* weights_Y2N(ind_tail)) );
Z1_tail(1,3) = max( abs(q2(ind_tail)) ./ (Lap2N(ind_tail) .* weights_Y2N(ind_tail)) );

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
D2G = zeros(3,3,3);
% the mu*psi term 
D2G(2,1,2) = 1;
D2G(2,2,1) = 1;
D2G(3,1,3) = 1;
D2G(3,3,1) = 1;

% norms of A (lazy and not very sharp estimate)
normA = opnorm_blocks_SKT_psi(Abar, nu, s-2, s-2, N, N);
normA(2:3,2:3) = normA(2:3,2:3) + cstemax*norm2x2(w, nu, s-2);

Z2 = reshape( normA' * reshape(D2G, [3, 9]), [3, 3, 3]);

% % Displaying the obtained estimates
% fprintf("\nZ2 =\n")
% disp(Z2)


%% Error terms coming from the fact that G involves DuF at the exact solution

% Remember that the error bounds in rX are in the strong norm (with n^k
% weights)
rlambda = rX(1);
ru1 = rX(2);
ru2 = rX(3);

beta = [rlambda+2*d11*ru1+d12*ru2, d12*ru1;
        d21*ru2, rlambda+d21*ru1+2*d22*ru2]; 
S = [2*a1*ru1+b1*ru2, b1*ru1;
     b2*ru2, b2*ru1+2*a2*ru2];

% The error terms (coming from DuF(X)-DuF(\bX))
errDFnoLap = zeros(3,3);
errDFLap = zeros(3,3);
if useintervals
    errDFnoLap = intval(errDFnoLap);
    errDFLap = intval(errDFLap);
end
errDFnoLap(2:3,2:3) = S;
errDFLap(2:3,2:3) = cstemin * beta;

% norm of A from k-2 to k (lazy and not very sharp estimate)
normLapA = opnorm_blocks_SKT_psi(Abar, nu, s, s-2, N, N);
normLapA(2:3,2:3) = normLapA(2:3,2:3) + 4*cstemax*norm2x2(w, nu, s);

% DFA
errAadjepsadj = ( errDFnoLap*normA + errDFLap*normLapA )';

% Updated bounds
weights_Y = getweights(N, nu, s-2);
norm_psi1 = max(abs(Psi1)./weights_Y);
norm_psi2 = max(abs(Psi2)./weights_Y);
Yerr = errAadjepsadj * [0; norm_psi1; norm_psi2];
Y = Y + Yerr;

Z1err = errAadjepsadj;
Z1 = Z1 + Z1err;

%% Checking that we have a contraction
rmin = polynomialsnegative(Y,Z1,Z2);


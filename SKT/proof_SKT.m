function [rmin, Abar, w, gamma] = proof_SKT(X, para, N, nsym, nu, s, Abar, w, gamma)

% Tries to prove the bifurcation point described by X =(lambda; usym; phi),
% i.e., to prove that the extended system has a zero near X.
% The last 3 arguments are optional, and are computed inside this function
% if not given as inputs.

%% Initialization

% Check whether we are using Intlab or not
useintervals = isa(X(1), 'intval');
if useintervals
    fprintf('\nRigorous validation with interval arithmetic\n')
    ipi = intval('pi');
else
    fprintf('\nPrevalidation without interval arithmetic\n')
    ipi = pi;
end

% check whether the solution is 1D or 2D
is1D = not(isfield(para,'Ly')); 

% get u_i and phi_i
[lambda, u1, u2, phi1, phi2] = extract_SKT(X, N, nsym, 'mat');

% This should already be true, but better be safe
para.d1 = lambda;
para.d2 = lambda;

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

% Larger index sets for computing the nonlinear terms exactly
N2 = 2*N-1;
X2N = projection_SKT(X, N2, N, nsym);

% total sizes
NN = prod(N);
NN2 = prod(N2);

% The vector of Fourier coefficients representing the constant function
e = zeros(N);
e(1,1) = 1;

% The approximate inverse A is defined by Abar, w and gamma
if nargin < 7
    % If Abar, w and gamma are not given, they are computed here
    DDif11 = d1*e + 2*d11*u1 + d12*u2;
    DDif12 = d12*u1;
    DDif21 = d21*u2;
    DDif22 = d2*e + d21*u1 + 2*d22*u2;

    MDDif = [convomat(DDif11), convomat(DDif12);
             convomat(DDif21), convomat(DDif22)];
    
    E = e(:);
    [~, indcomplN] = getindsym(N, nsym);

    W_dot1 = MDDif \ [E;zeros(NN,1)];
    W_dot1(indcomplN) = 0;
    W11 = W_dot1(1:NN);
    W21 = W_dot1(NN+(1:NN));

    W_dot2 = MDDif \ [zeros(NN,1);E];
    W_dot2(indcomplN) = 0;
    W12 = W_dot2(1:NN);
    W22 = W_dot2(NN+(1:NN));

    w = {reshape(W11,N), reshape(W12,N);
         reshape(W21,N), reshape(W22,N)};

    beta = {2*d11*phi1+d12*phi2, d12*phi1;
            d21*phi2, d21*phi1+2*d22*phi2};
    
    gamma = myprod2x2(w, myprod2x2(beta,w));
    gamma = {-gamma{1,1},-gamma{1,2};
             -gamma{2,1},-gamma{2,2}};

    [Fext2N, DFext2N] = F_DF_ext_SKT(X2N, para, N2, nsym);
    Abar = inv(DFext2N);

    % The indices to extract only the part of size N from the whole 2N thing
    ind_Nin2N = getindwithin(N, N2);
    [ind_NsyminN2sym, ~, ~, N2sym] = getindsymwithin(N, N2, nsym);
    ind_XNinX2N = [1; 1+(ind_NsyminN2sym); 1+N2sym+(ind_NsyminN2sym); 1+2*N2sym+ind_Nin2N; 1+2*N2sym+NN2+ind_Nin2N];
    
    Abar = Abar(ind_XNinX2N, ind_XNinX2N);
else
    [~, ~, N2sym] = getindsym(N2, nsym);
    Fext2N = F_DF_ext_SKT(X2N, para, N2, nsym);
end

%% Y
N3 = 3*N-2;
AF = applyA_SKT(Abar, w, gamma, Fext2N, para, N, N2, N3, nsym);
Y = norm_blocks_SKT(AF, nu, s, N3, nsym);

% Displaying the obtained estimates
fprintf("\nY =\n")
disp(Y)

%% Z1
%%% Finite part
N4 = 4*N-3;
X3N = projection_SKT(X, N3, N, nsym);
% Square projection of DFext (the largest one we need)
[~,DFext3N] = F_DF_ext_SKT(X3N, para, N3, nsym);

ind_N2inN3 = getindwithin(N2, N3);
[ind_N2syminN3sym, ~, ~, N3sym] = getindsymwithin(N2, N3, nsym);
NN3 = prod(N3);
ind_X2NinX3N = [1; 1+(ind_N2syminN3sym); 1+N3sym+(ind_N2syminN3sym); 1+2*N3sym+ind_N2inN3; 1+2*N3sym+NN3+ind_N2inN3];
% Getting only the required columns of DFext 
DFext3N = DFext3N(:,ind_X2NinX3N); 
% The finite part of A*DFext
ADF = applyA_SKT(Abar, w, gamma, DFext3N, para, N, N3, N4, nsym);

% Constructing the identity matrix (non square)
[ind_N2syminN4sym, ~, ~, N4sym] = getindsymwithin(N2, N4, nsym);
ind_N2inN4 = getindwithin(N2, N4);
NN4 = prod(N4);
N2full = 1 + 2*N2sym + 2*NN2;
N4full = 1 + 2*N4sym + 2*NN4;
I_4N2N = zeros(N4full, N2full);
I_4N2N(1) = 1;
I_4N2N(1+ind_N2syminN4sym,1+(1:N2sym)) = eye(N2sym);
I_4N2N(1+N4sym+ind_N2syminN4sym,1+N2sym+(1:N2sym)) = eye(N2sym);
I_4N2N(1+2*N4sym+ind_N2inN4,1+2*N2sym+(1:NN2)) = eye(NN2);
I_4N2N(1+2*N4sym+NN4+ind_N2inN4,1+2*N2sym+NN2+(1:NN2)) = eye(NN2);

% Finite part of Z1
B = I_4N2N - ADF;
Z1_finite = opnorm_blocks_SKT(B, nu, s, s, N4, N2, nsym);

% % Displaying the obtained estimates
% fprintf("\nZ1_finite =\n")
% disp(Z1_finite)

%%% tail part 
eN2 = zeros(N2);
eN2(1,1) = 1;
alpha = {d1*e+2*d11*u1+d12*u2, d12*u1;
         d21*u2, d2*e+d21*u1+2*d22*u2}; 
beta = {2*d11*phi1+d12*phi2, d12*phi1;
        d21*phi2, d21*phi1+2*d22*phi2};
R = {r1*e-2*a1*u1-b1*u2, -b1*u1;
     -b2*u2, r2*e-b2*u1-2*a2*u2}; 
S = {-2*a1*phi1-b1*phi2, -b1*phi1;
     -b2*phi2, -b2*phi1-2*a2*phi2};
delta = myprod2x2(w,alpha,N2);
delta = {eN2-delta{1,1},-delta{1,2};
         -delta{2,1},eN2-delta{2,2}};

eta1 = myprod2x2(gamma, alpha, N2);
eta2 = myprod2x2(w, beta, N2);
eta = {eta1{1,1}+eta2{1,1}, eta1{1,2}+eta2{1,2};
       eta1{2,1}+eta2{2,1}, eta1{2,2}+eta2{2,2}};

normw = norm2x2(w, nu, s);
normgamma = norm2x2(gamma, nu, s);
normR = norm2x2(R,  nu, s);
normS = norm2x2(S, nu, s);
normdelta = norm2x2(delta, nu, s);
normeta = norm2x2(eta, nu, s);

Z1_tail = zeros(5,5);
if useintervals
    Z1_tail = intval(Z1_tail);
end
if is1D
    lambdaN = (ipi/L1)^2*N(1)^2;
else
    lambdaN = min((ipi/L1)^2*N(1)^2, (ipi/L2)^2*N(2)^2);
end
Z1_tail(2:3,2:3) = normdelta + 1/lambdaN * normw * normR;
Z1_tail(4:5,2:3) = normeta + 1/lambdaN * (normgamma*normR + normw*normS);
Z1_tail(4:5,4:5) = Z1_tail(2:3,2:3);

% % Displaying the obtained estimates
% fprintf("\nZ1_tail =\n")
% disp(Z1_tail)
   
%%% full Z1
Z1 = max(Z1_finite,Z1_tail);

% Displaying the obtained estimates
fprintf("\nZ1 =\n")
disp(Z1)
fprintf("\nSpectral radius of Z_1 : %f \n",eigs(i2f(Z1), 1)) % If this is >= 1, the proof will fail (it may also fail if this is too close to 1)

   
%% Z2
% /!\ Important remark: For Z2 we actually do something slightly different
% from what is written in the paper. That is, instead of using estimates of
% the form ||A*D^2F||_{X->X} <= ||A||_{Y->X} * ||D^2F||_{X->Y}, we split
% D^2F between bounded and unbounded terms. For the bounded ones, we simply
% use ||A*D^2F||_{X->X} <= ||A||_{X->X} * ||D^2F||_{X->X}. For the
% unbounded ones, we do instead
% ||A*D^2F||_{X->X} <= ||A*\Delta||_{X->X} * ||\Delta^{-1}*D^2F||_{X->X}.

D2FLap = zeros(5,5,5);
D2FnoLap = zeros(5,5,5);
if useintervals
    D2FLap = intval(D2FLap);
    D2FnoLap = intval(D2FnoLap);
end

% DlambdaF
for i = 2:5
    D2FLap(i,1,i) = 1;
end

% DuF (sym)
D2FLap(2:3,2:3,1) = [1 0;
                     0 1];
D2FLap(2:3,2:3,2) = [2*d11  d12;
                     0      d21];
D2FLap(2:3,2:3,3) = [d12  0;
                     d21  2*d22];
D2FnoLap(2:3,2:3,2) = [2*a1  b1;
                       0     b2];
D2FnoLap(2:3,2:3,3) = [b1  0;
                       b2  2*a2];

% DuuFphi
D2FLap(4:5,2:3,4) = [2*d11  d12;
                     0      d21];
D2FLap(4:5,2:3,5) = [d12  0;
                     d21  2*d22];
D2FnoLap(4:5,2:3,4) = [2*a1  b1;
                       0     b2];
D2FnoLap(4:5,2:3,5) = [b1  0;
                       b2  2*a2];

% DuF (nosym)
D2FLap(4:5,4:5,1) = [1 0;
                     0 1];
D2FLap(4:5,4:5,2) = [2*d11  d12;
                     0      d21];
D2FLap(4:5,4:5,3) = [d12  0;
                     d21  2*d22];
D2FnoLap(4:5,4:5,2) = [2*a1  b1;
                       0     b2];
D2FnoLap(4:5,4:5,3) = [b1  0;
                       b2  2*a2];

% norms of A and ALap (lazy and not very sharp estimate)
normA = opnorm_blocks_SKT(Abar, nu, s, s, N, N, nsym);
if is1D
    cste = 1/(ipi/L1)^2;
else
    cste = 1/(ipi/max(L1,L2))^2;
end
normA(2:3,2:3) = normA(2:3,2:3) + cste*normw;
normA(4:5,2:3) = normA(4:5,2:3) + cste*normgamma;
normA(4:5,4:5) = normA(4:5,4:5) + cste*normw;

if is1D
    Lapu = (ipi/L1)^2 * ((0:N(1)-1)').^2;
else
    [N2, N1] = meshgrid(0:N(2)-1,0:N(1)-1);
    Lapu = (ipi/L1)^2*N1.^2 + (ipi/L2)^2*N2.^2;
end
indsymN = getindsym(N, nsym);
LapX = [0;Lapu(indsymN);Lapu(indsymN);Lapu(:);Lapu(:)]';
normALap = opnorm_blocks_SKT(Abar.*LapX, nu, s, s, N, N, nsym);
normALap(2:3,2:3) = normALap(2:3,2:3) + normw;
normALap(4:5,2:3) = normALap(4:5,2:3) + normgamma;
normALap(4:5,4:5) = normALap(4:5,4:5) + normw;

Z2Lap = reshape( normALap * reshape(D2FLap, [5, 25]), [5, 5, 5]);
Z2noLap = reshape( normA * reshape(D2FnoLap, [5, 25]), [5, 5, 5]);
Z2 = Z2Lap + Z2noLap;

% % Displaying the obtained estimates
% fprintf("\nZ2 =\n")
% disp(Z2)

%% Checking that we have a contraction
rmin = polynomialsnegative(Y,Z1,Z2);


function [rmin, Abar] = proof_OK(X, para, N, nsym, nu, s, rstar, Abar)

% Tries to prove the bifurcation point described by X=(lambda, usym, phi).
% The last argument is optional, and is computed inside this function if 
% not given as input.

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

% get lambda, u and phi
[lambda, u, phi] = extract_OK(X, N, nsym, 'mat');

% This should already be true, but better be safe
para.lambda = lambda;

% All the parameters of the system
sigma = para.sigma;
L1 = para.Lx;
if not(is1D)
    L2 = para.Ly;
end

% The approximate inverse A is defined by Abar
if nargin < 8
    % If Abar is not given, it is computed here
    [~, DFext] = F_DF_ext_OK(X, para, N, nsym);
    Abar = inv(DFext);
end


%% Y
N3 = 3*N-2;
X3N = projection_OK(X, N3, N, nsym);
Fext3N = F_DF_ext_OK(X3N, para, N3, nsym);
AF = applyA_OK(Abar, Fext3N, para, N, N3, nsym);
Y = norm_blocks_OK(AF, nu, s, N3, nsym);

% Displaying the obtained estimates
fprintf("\nY =\n")
disp(Y)

%% Z1
%%% Finite part
N5 = 5*N-4;
X5N = projection_OK(X, N5, N, nsym);
% Square projection of DFext (the largest one we need)
[~,DFext5N] = F_DF_ext_OK(X5N, para, N5, nsym);

ind_N3inN5 = getindwithin(N3, N5);
[ind_N3syminN5sym, ~, ~, N5sym] = getindsymwithin(N3, N5, nsym);
ind_X3NinX5N = [1; 1+(ind_N3syminN5sym); 1+N5sym+ind_N3inN5];
% Getting only the required columns of DFext
DFext5N = DFext5N(:,ind_X3NinX5N); 
% The finite part of A*DFext
ADF = applyA_OK(Abar, DFext5N, para, N, N5, nsym);

% Constructing the identity matrix (non square)
NN3 = prod(N3);
NN5 = prod(N5);
[~, ~, N3sym] = getindsym(N3, nsym);
N3full = 1 + N3sym + NN3;
N5full = 1 + N5sym + NN5;
I_5N3N = zeros(N5full, N3full);
I_5N3N(1) = 1;
I_5N3N(1+ind_N3syminN5sym,1+(1:N3sym)) = eye(N3sym);
I_5N3N(1+N5sym+ind_N3inN5,1+N3sym+(1:NN3)) = eye(NN3);

% Finite part of Z1
B = I_5N3N - ADF;
Z1_finite = opnorm_blocks_OK(B, nu, s, s, N5, N3, nsym);

% % Displaying the obtained estimates
% fprintf("\nZ1_finite =\n")
% disp(Z1_finite)

%%% tail part
N2 = 2*N-1;
eN2 = zeros(N2);
eN2(1,1) = 1;

weightsN2 = getweights(N2, nu, s);

Mu = convomat(u, N2, N);
norm_1m3u2 = weightsN2' * abs(eN2(:) - 3*Mu*u(:));
norm_uphi = weightsN2' * abs(Mu*phi(:));

Z1_tail = zeros(3,3);
if useintervals
    Z1_tail = intval(Z1_tail);
end
if is1D
    lambdaN = (ipi/L1)^2*N(1)^2;
    lambda3N = (ipi/L1)^2*N3(1)^2;
else
    lambdaN = min((ipi/L1)^2*N(1)^2, (ipi/L2)^2*N(2)^2);
    lambda3N = min((ipi/L1)^2*N3(1)^2, (ipi/L2)^2*N3(2)^2);
end
Z1_tail(2,2) = abs(lambda) * (norm_1m3u2/lambdaN +abs(sigma)/lambda3N^2);
Z1_tail(3,2) = 6*abs(lambda) * norm_uphi / lambdaN;
Z1_tail(3,3) = Z1_tail(2,2);

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
% than what is written in the paper. That is, instead of using estimates of
% the form ||A*D^2F||_{X->X} <= ||A||_{Y->X} * ||D^2F||_{X->Y}, we split
% D^2F between bounded and unbounded terms. For the bounded ones, we simply
% use ||A*D^2F||_{X->X} <= ||A||_{X->X} * ||D^2F||_{X->X}. For the
% unbounded ones, we do instead
% ||A*D^2F||_{X->X} <= ||A*\Delta||_{X->X} * ||\Delta^{-1}*D^2F||_{X->X}.

weightsN = getweights(N, nu, s);
norm_u = weightsN' * abs(u(:));
norm_phi = weightsN' * abs(phi(:));

norm_1m3u2_star = norm_1m3u2 + 3*(2*norm_u*rstar+rstar^2);
norm_lambdau_star = (abs(lambda)+rstar) * (norm_u+rstar);
norm_uphi_star = (norm_u+rstar) * (norm_phi+rstar);
norm_lambdaphi_star = (abs(lambda)+rstar) * (norm_phi+rstar);

D2FLap = zeros(3,3,3);
D2FnoLap = zeros(3,3,3);
if useintervals
    D2FLap = intval(D2FLap);
    D2FnoLap = intval(D2FnoLap);
end

% DlambdaF
D2FLap(2,1,2) = norm_1m3u2_star;
D2FnoLap(2,1,2) = abs(sigma);

% DlambdauFphi
D2FLap(3,1,2) = 6*norm_uphi_star;
D2FLap(3,1,3) = norm_1m3u2_star;
D2FnoLap(3,1,3) = abs(sigma);

% DuF (sym)
D2FLap(2,2,1) = norm_1m3u2_star;
D2FLap(2,2,2) = 6*norm_lambdau_star;
D2FnoLap(2,2,1) = abs(sigma);

% DuuFphi
D2FLap(3,2,1) = 6*norm_uphi_star;
D2FLap(3,2,2) = 6*norm_lambdaphi_star;
D2FLap(3,2,3) = 6*norm_lambdau_star;

% DuF (nosym)
D2FLap(3,3,1) = norm_1m3u2_star;
D2FLap(3,3,2) = 6*norm_lambdau_star;
D2FnoLap(3,3,1) = abs(sigma);

% norms of A and ALap
normA = opnorm_blocks_OK(Abar, nu, s, s, N, N, nsym);
normA(2,2) = max(normA(2,2), 1/lambdaN^2);
normA(3,3) = max(normA(3,3), 1/lambdaN^2);

if is1D
    Lapu = (ipi/L1)^2 * ((0:N(1)-1)').^2;
else
    [N2, N1] = meshgrid(0:N(2)-1,0:N(1)-1);
    Lapu = (ipi/L1)^2*N1.^2 + (ipi/L2)^2*N2.^2;
end
indsymN = getindsym(N, nsym);
LapX = [0;Lapu(indsymN);Lapu(:)]';
normALap = opnorm_blocks_OK(Abar.*LapX, nu, s, s, N, N, nsym);
normALap(2,2) = max(normALap(2,2), 1/lambdaN);
normALap(3,3) = max(normALap(3,3), 1/lambdaN);

Z2Lap = reshape( normALap * reshape(D2FLap, [3, 9]), [3, 3, 3]);
Z2noLap = reshape( normA * reshape(D2FnoLap, [3, 9]), [3, 3, 3]);
Z2 = Z2Lap + Z2noLap;

% % Displaying the obtained estimates
% fprintf("\nZ2 =\n")
% disp(Z2)

%% Checking that we have a contraction
rmin = polynomialsnegative(Y,Z1,Z2);
if any(rstar<rmin)
    rmin = NaN;
end


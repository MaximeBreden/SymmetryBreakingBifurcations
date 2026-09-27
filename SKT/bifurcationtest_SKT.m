function [psiDuuF, rpsiDuuF] = bifurcationtest_SKT(X, Z, para, N, nsym, nu, s, rX, rZ)

% Computes psi^* [D_{uu}F(lambda,u)(phi,phi)], which is zero iff the bifurcation
% is a pitchfork. X = (lambda; usym; phi) and Z = (mu; psi).
% If error bounds rX and rZ for X and Z are provided, then also computes an
% error bound for psi^* [D_{uu}F(X)(phi,phi)].

%% Initialization
% Check whether we are using Intlab and intervals or not
if isa(X(1), 'intval')
    ipi = intval('pi');
else
    ipi = pi;
end

% check whether the solution is 1D or 2D
is1D = not(isfield(para,'Ly')); 

[~, ~, ~, phi1, phi2] = extract_SKT(X, N, nsym, 'mat');
[~, Psi1, Psi2] = extract_SKT_psi(Z, N, 'vect');

% the required parameters of the system
d11 = para.d11;
d12 = para.d12;
d21 = para.d21;
d22 = para.d22;
a1 = para.a1;
b1 = para.b1;
a2 = para.a2;
b2 = para.b2;
Lx = para.Lx;
if not(is1D)
    Ly = para.Ly;
end

% Building the Laplacian operator
N2 = 2*N-1;
if is1D
    Lap = -(ipi/Lx)^2 * ((0:N2(1)-1)').^2;
else
    [N2y, N2x] = meshgrid(0:N2(2)-1,0:N2(1)-1);
    Lap = -(ipi/Lx)^2*N2x.^2 - (ipi/Ly)^2*N2y.^2;
end

%% Computing psi^*(DuuF(lambda,u)(phi,phi))
phi1phi1 = convo(phi1, phi1, N2);
phi1phi2 = convo(phi1, phi2, N2);
phi2phi2 = convo(phi2, phi2, N2);
LapPhi1Phi1 = Lap(:) .* phi1phi1(:);
LapPhi1Phi2 = Lap(:) .* phi1phi2(:);
LapPhi2Phi2 = Lap(:) .* phi2phi2(:);

DuuF1 = 2 * ( d11*LapPhi1Phi1 + d12*LapPhi1Phi2 - a1*phi1phi1(:) - b1*phi1phi2(:) );
DuuF2 = 2 * ( d21*LapPhi1Phi2 + d22*LapPhi2Phi2 - b2*phi1phi2(:) - a2*phi2phi2(:) );

indNin2N = getindwithin(N, N2);
psiDuuF = sum( Psi1.*DuuF1(indNin2N) + Psi2.*DuuF2(indNin2N) );

%% Adding the error term
if nargin > 5
    if s < 0
        error("The error estimates that are implemented assume that \ell^1_{s} is a Banach algebra, and hence that s>=0.")
    end

    rphi1 = rX(4);
    rphi2 = rX(5);
    rpsi1 = rZ(2);
    rpsi2 = rZ(3);
    
    % error bounds for the phiiphij terms in X = \ell^1_{s}
    weights_X = getweights(N, nu, s);
    normphi1_X = weights_X' * abs(phi1(:));
    normphi2_X = weights_X' * abs(phi2(:));

    rphi1phi1_X = 2*normphi1_X*rphi1 + rphi1^2;
    rphi1phi2_X = normphi1_X*rphi2 + normphi2_X*rphi1 + rphi1*rphi2;
    rphi2phi2_X = 2*normphi2_X*rphi2 + rphi2^2;

    % error bounds for the phiiphij terms in Y = \ell^1_{s-2}
    % (||.||_Y <= ||.||_X, hence the error bounds rphi, which are in the X norm, can still be used)
    weights_Y = getweights(N, nu, s-2);
    normphi1_Y = weights_Y' * abs(phi1(:));
    normphi2_Y = weights_Y' * abs(phi2(:));

    rphi1phi1_Y = 2*normphi1_Y*rphi1 + rphi1^2;
    rphi1phi2_Y = normphi1_Y*rphi2 + normphi2_Y*rphi1 + rphi1*rphi2;
    rphi2phi2_Y = 2*normphi2_Y*rphi2 + rphi2^2;

    % error bounds for DuuF, in Y = \ell^1_{k-2}
    if is1D
        cste = (ipi/Lx)^2;
    else
        cste = (ipi/min(Lx,Ly))^2;
    end
    rDuuF1 = 2*( cste * (d11*rphi1phi1_X+d12*rphi1phi2_X) + a1*rphi1phi1_Y+b1*rphi1phi2_Y );
    rDuuF2 = 2*( cste * (d21*rphi1phi2_X+d22*rphi2phi2_X) + b2*rphi1phi2_Y+a2*rphi2phi2_Y );

    % error bound for psi^*(DuuF(lambda,u)(phi,phi))
    normpsi1 = max( abs(Psi1) ./ weights_Y );
    normpsi2 = max( abs(Psi2) ./ weights_Y );
    weights_Y_2N = getweights(N2, nu, s-2)';
    normDuuF1 = weights_Y_2N * abs(DuuF1);
    normDuuF2 = weights_Y_2N * abs(DuuF2);

    rpsiDuuF = normpsi1*rDuuF1+normpsi2*rDuuF2 + rpsi1*normDuuF1+rpsi2*normDuuF2 + rpsi1*rDuuF1+rpsi2*rDuuF2;
end
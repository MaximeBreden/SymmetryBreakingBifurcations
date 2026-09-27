function [psiDuuF, rpsiDuuF] = bifurcationtest_OK(X, Z, para, N, nsym, nu, s, rX, rZ)

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

[lambda, u, phi] = extract_OK(X, N, nsym, 'mat');
[~, Psi] = extract_OK_psi(Z, N);

% the required parameters of the system
Lx = para.Lx;
if not(is1D)
    Ly = para.Ly;
end

%
N2 = 2*N-1;
N3 = 3*N-2;

% Building the Laplacian operator
if is1D
    Lap = -(ipi/Lx)^2 * ((0:N3(1)-1)').^2;
else
    [N3y, N3x] = meshgrid(0:N3(2)-1,0:N3(1)-1);
    Lap = -(ipi/Lx)^2*N3x.^2 - (ipi/Ly)^2*N3y.^2;
end

%% Computing psi^*(DuuF(lambda,u)(phi,phi)) = 6*lambda * psi^*(Delta(u*phi^2))
phiphi = convo(phi, phi, N2);
uphiphi = convo(u, phiphi, N3);
Lapuphiphi = Lap(:) .* uphiphi(:);
indNin3N = getindwithin(N, N3);
psiLapuphiphi = sum( Psi .* Lapuphiphi(indNin3N) );
psiDuuF = 6 * lambda * psiLapuphiphi;

%% Adding the error term
if nargin > 5
    if s < 2
        error("The upcoming error estimate assumes that \ell^1_{k-2} is a Banach algebra, and hence that k>=2.")
    end

    rlambda = rX(1);
    ru = rX(2);
    rphi = rX(3);
    rpsi = rZ(2);

    normu = sum( abs(u(:)) .* getweights(N, nu, s-2) );
    normphi = sum( abs(phi(:)) .* getweights(N, nu, s-2) );
    normphiphi = sum( abs(phiphi(:)) .* getweights(N2, nu, s-2) );

    % error bound for Delta(u*phi^2), in \ell^1_{s-4}
    if is1D
        cste = (ipi/Lx)^2;
    else
        cste = (ipi/min(Lx,Ly))^2;
    end
    rLapuphiphi = cste * ( ru*(normphiphi+2*rphi*normphi+rphi^2) + normu*(2*rphi*normphi+rphi^2) );

    normpsi = max( abs(Psi) ./ getweights(N, nu, s-4));
    normLapuphiphi = sum( abs(Lapuphiphi) .* getweights(N3, nu, s-4));

    % error bound for psi^*(Delta(u*phi^2))
    rpsiLapuphiphi = normpsi*rLapuphiphi + rpsi*(normLapuphiphi+rLapuphiphi);

    rpsiDuuF = 6 * ( abs(lambda)*rpsiLapuphiphi + rlambda*abs(psiLapuphiphi) + rlambda*rpsiLapuphiphi );
end
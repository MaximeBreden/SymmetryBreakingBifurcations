function Y = applyAadj_OK_psi(Abar, Z, para, NA, NZ)

% Computes Y = A^* * Z, without truncating
% Z = [mu; psi], and Z can have multiple columns.

%% Initialization

% Check whether we are using Intlab and intervals or not
if isa(Z(1), 'intval')
    useintervals = true;
    ipi = intval('pi');
else
    useintervals = false;
    ipi = pi;
end

% check whether the solution is 1D or 2D
is1D = not(isfield(para,'Ly'));

NNA = prod(NA);
NNZ = prod(NZ);

% Inverse Laplacian squared
if is1D
    InvLap2 = 1 ./ (-(ipi/para.Lx)^2 * ((0:NZ(1)-1)').^2).^2;
else
    [N2, N1] = meshgrid(0:NZ(2)-1,0:NZ(1)-1);
    InvLap2 = 1./(-(ipi/para.Lx)^2*N1.^2 - (ipi/para.Ly)^2*N2.^2).^2;
    InvLap2 = InvLap2(:);
end
InvLap2(1,1) = 0;

%% Lots of index sets, to be able to extract the various parts when computing A*Z

[ind_psiAinpsiZ, ind_notpsiAinpsiZ] = getindwithin(NA, NZ);
ind_AinZ = [1; 1+(ind_psiAinpsiZ)]; 
ind_psiAinA = 1+(1:NNA);

%% precomputation of the tails (the \Delta^{-2} acting outside of Abar)
[~, Psi] = extract_OK_psi(Z, NZ, 'vect');
InvLap2psi_notAinZ = InvLap2(ind_notpsiAinpsiZ) .* Psi(ind_notpsiAinpsiZ,:);

%% Y = A^* * Z
Abaradj = Abar';

Ymu = Abaradj(1,:) * Z(ind_AinZ,:);

Ypsi = zeros(NNZ,size(Z,2));
if useintervals
    Ypsi = intval(Ypsi);
end
Ypsi(ind_psiAinpsiZ,:) = Abaradj(ind_psiAinA,:)*Z(ind_AinZ,:);
Ypsi(ind_notpsiAinpsiZ,:) = -InvLap2psi_notAinZ;

Y = [Ymu; 
     Ypsi];
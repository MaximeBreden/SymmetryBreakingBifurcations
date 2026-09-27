function Y = applyAadj_SKT_psi(Abar, w, Z, para, NA, NZ, NY)

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
NNY = prod(NY);

if is1D
    InvLap = 1 ./ (-(ipi/para.Lx)^2 * ((0:NY(1)-1)').^2);
else
    [N2, N1] = meshgrid(0:NY(2)-1,0:NY(1)-1);
    InvLap = 1./(-(ipi/para.Lx)^2*N1.^2 - (ipi/para.Ly)^2*N2.^2);
    InvLap = InvLap(:);
end
InvLap(1,1) = 0;

% /!\ We incorporate the inverse Laplacian in the Mwij matrices
Mw11 = InvLap .* convomat(w{1,1}, NY, NZ);
Mw12 = InvLap .* convomat(w{1,2}, NY, NZ);
Mw21 = InvLap .* convomat(w{2,1}, NY, NZ);
Mw22 = InvLap .* convomat(w{2,2}, NY, NZ);

%% Lots of index sets, to be able to extract the various parts when computing A*Z

[ind_uAinuZ, ind_notuAinuZ] = getindwithin(NA, NZ);

ind_AinZ = [1; 1+(ind_uAinuZ); 1+NNZ+(ind_uAinuZ)]; 
ind_psiAinA = 1+(1:2*NNA);

[ind_uAinuY, ind_notuAinuY] = getindwithin(NA, NY);
ind12_uAinuY = [ind_uAinuY; NNY+ind_uAinuY];
ind12_notuAinuY = [ind_notuAinuY; NNY+ind_notuAinuY];


%% Y = A^* * Z
Abaradj = Abar';

% Ymu
Ymu = Abaradj(1,:) * Z(ind_AinZ,:);

% Ypsi
psi1 = Z(1+(1:NNZ),:);
psi2 = Z(1+NNZ+(1:NNZ),:);
Ypsi = zeros(2*NNY,size(Z,2));
if useintervals
    Ypsi = intval(Ypsi);
end
Ypsi(ind12_uAinuY,:) = Abaradj(ind_psiAinA,:)*Z(ind_AinZ,:) + ...
                       [Mw11(ind_uAinuY,ind_notuAinuZ)*psi1(ind_notuAinuZ,:) + Mw21(ind_uAinuY,ind_notuAinuZ)*psi2(ind_notuAinuZ,:);
                        Mw12(ind_uAinuY,ind_notuAinuZ)*psi1(ind_notuAinuZ,:) + Mw22(ind_uAinuY,ind_notuAinuZ)*psi2(ind_notuAinuZ,:)];
Ypsi(ind12_notuAinuY,:) = [Mw11(ind_notuAinuY,:)*psi1 + Mw21(ind_notuAinuY,:)*psi2;
                           Mw12(ind_notuAinuY,:)*psi1 + Mw22(ind_notuAinuY,:)*psi2];

Y = [Ymu; 
     Ypsi];
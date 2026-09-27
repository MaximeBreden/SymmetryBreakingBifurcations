function Y = applyA_SKT(Abar, w, gamma, X, para, NA, NX, NY, nsym)

% Computes Y = A*X, without truncating.
% X = [lambda; usym; phi], and X can have multiple columns.

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

NNA = prod(NA);
NNX = prod(NX);
NNY = prod(NY);

% Inverse Laplacian
if is1D
    InvLap = 1 ./ (-(ipi/para.Lx)^2 * ((0:NX(1)-1)').^2);
else
    [N2, N1] = meshgrid(0:NX(2)-1,0:NX(1)-1);
    InvLap = 1./(-(ipi/para.Lx)^2*N1.^2 - (ipi/para.Ly)^2*N2.^2);
    InvLap = InvLap(:);
end
InvLap(1,1) = 0;

Mw11 = convomat(w{1,1}, NY, NX);
Mw12 = convomat(w{1,2}, NY, NX);
Mw21 = convomat(w{2,1}, NY, NX);
Mw22 = convomat(w{2,2}, NY, NX);

Mgamma11 = convomat(gamma{1,1}, NY, NX);
Mgamma12 = convomat(gamma{1,2}, NY, NX);
Mgamma21 = convomat(gamma{2,1}, NY, NX);
Mgamma22 = convomat(gamma{2,2}, NY, NX);

%% Lots of index sets, to be able to extract the various parts when computing A*X

[ind_usymXinuX, ~, NsymX] = getindsym(NX, nsym); 
ind_usymXinX = 1+(1:2*NsymX);
ind_phiinX = 1+2*NsymX+(1:2*NNX);

[ind_usymYinuY, ~, NsymY] = getindsym(NY, nsym);

[~, ~, NsymA] = getindsym(NA, nsym);
ind_usymAinA = 1+(1:2*NsymA);

[ind_usymAinusymX, ind_notusymAinusymX, ind_notuAinuX_sym] = getindsymwithin(NA, NX, nsym);
[ind_uAinuX, ind_notuAinuX] = getindwithin(NA, NX);

ind_AinX = [1; 1+(ind_usymAinusymX); 1+NsymX+(ind_usymAinusymX); 1+2*NsymX+(ind_uAinuX); 1+2*NsymX+NNX+(ind_uAinuX)]; 
ind_phiAinA = 1+2*NsymA+(1:2*NNA);

[ind_uAinuY, ind_notuAinuY] = getindwithin(NA, NY);
ind12_uAinuY = [ind_uAinuY; NNY+ind_uAinuY];
[ind_usymAinusymY, ind_notusymAinusymY, ind_notuAinuY_sym] = getindsymwithin(NA, NY, nsym);
ind12_usymAinusymY = [ind_usymAinusymY; NsymY+ind_usymAinusymY];
ind12_notusymAinusymY = [ind_notusymAinusymY; NsymY+ind_notusymAinusymY];
ind_usymAinuY = ind_usymYinuY(ind_usymAinusymY);
ind12_notuAinuY = [ind_notuAinuY; NNY+ind_notuAinuY];

%% Precomputation
usym = X(ind_usymXinX,:);
phi = X(ind_phiinX,:);

InvLapu1sym = InvLap(ind_usymXinuX) .* usym(1:NsymX,:);
InvLapu2sym = InvLap(ind_usymXinuX) .* usym(NsymX+(1:NsymX),:);
InvLapphi1 = InvLap .* phi(1:NNX,:);
InvLapphi2 = InvLap .* phi(NNX+(1:NNX),:);

% clear usym phi

%% Y = A*X

% lambda part
Ylambda = Abar(1,:) * X(ind_AinX,:);

% usym part
Yusym = zeros(2*NsymY,size(X,2));
if useintervals
    Yusym = intval(Yusym);
end
Yusym(ind12_usymAinusymY,:) = Abar(ind_usymAinA,:)*X(ind_AinX,:) + ...
                              [Mw11(ind_usymAinuY,ind_notuAinuX_sym)*InvLapu1sym(ind_notusymAinusymX,:) + Mw12(ind_usymAinuY,ind_notuAinuX_sym)*InvLapu2sym(ind_notusymAinusymX,:);
                               Mw21(ind_usymAinuY,ind_notuAinuX_sym)*InvLapu1sym(ind_notusymAinusymX,:) + Mw22(ind_usymAinuY,ind_notuAinuX_sym)*InvLapu2sym(ind_notusymAinusymX,:)];
Yusym(ind12_notusymAinusymY,:) = [Mw11(ind_notuAinuY_sym,ind_usymXinuX)*InvLapu1sym + Mw12(ind_notuAinuY_sym,ind_usymXinuX)*InvLapu2sym;
                                  Mw21(ind_notuAinuY_sym,ind_usymXinuX)*InvLapu1sym + Mw22(ind_notuAinuY_sym,ind_usymXinuX)*InvLapu2sym];

% phi part
Yphi = zeros(2*NNY,size(X,2));
if useintervals
    Yphi = intval(Yphi);
end
Yphi(ind12_uAinuY,:) = Abar(ind_phiAinA,:)*X(ind_AinX,:) + ...
                       [Mgamma11(ind_uAinuY,ind_notuAinuX_sym)*InvLapu1sym(ind_notusymAinusymX,:) + Mgamma12(ind_uAinuY,ind_notuAinuX_sym)*InvLapu2sym(ind_notusymAinusymX,:);
                        Mgamma21(ind_uAinuY,ind_notuAinuX_sym)*InvLapu1sym(ind_notusymAinusymX,:) + Mgamma22(ind_uAinuY,ind_notuAinuX_sym)*InvLapu2sym(ind_notusymAinusymX,:)] + ...
                       [Mw11(ind_uAinuY,ind_notuAinuX)*InvLapphi1(ind_notuAinuX,:) + Mw12(ind_uAinuY,ind_notuAinuX)*InvLapphi2(ind_notuAinuX,:);
                        Mw21(ind_uAinuY,ind_notuAinuX)*InvLapphi1(ind_notuAinuX,:) + Mw22(ind_uAinuY,ind_notuAinuX)*InvLapphi2(ind_notuAinuX,:)];
Yphi(ind12_notuAinuY,:) = [Mgamma11(ind_notuAinuY,ind_usymXinuX)*InvLapu1sym + Mgamma12(ind_notuAinuY,ind_usymXinuX)*InvLapu2sym;
                           Mgamma21(ind_notuAinuY,ind_usymXinuX)*InvLapu1sym + Mgamma22(ind_notuAinuY,ind_usymXinuX)*InvLapu2sym] + ...
                          [Mw11(ind_notuAinuY,:)*InvLapphi1 + Mw12(ind_notuAinuY,:)*InvLapphi2;
                           Mw21(ind_notuAinuY,:)*InvLapphi1 + Mw22(ind_notuAinuY,:)*InvLapphi2];

Y = [Ylambda; 
     Yusym; 
     Yphi];
function Y = applyA_OK(Abar, X, para, NA, NX, nsym)

% Computes Y = A*X, without truncating.
% X = [lambda; usym; phi], and X can have multiple columns.

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

% Inverse Laplacian squared
if is1D
    InvLap2 = 1 ./ (-(ipi/para.Lx)^2 * ((0:NX(1)-1)').^2).^2;
else
    [N2, N1] = meshgrid(0:NX(2)-1,0:NX(1)-1);
    InvLap2 = 1./(-(ipi/para.Lx)^2*N1.^2 - (ipi/para.Ly)^2*N2.^2).^2;
    InvLap2 = InvLap2(:);
end
InvLap2(1,1) = 0;

%% Lots of index sets, to be able to extract the various parts when computing A*X

[ind_usymXinuX, ~, NsymX] = getindsym(NX, nsym); 

[~, ~, NsymA] = getindsym(NA, nsym);
ind_usymAinA = 1+(1:NsymA);

[ind_usymAinusymX, ind_notusymAinusymX] = getindsymwithin(NA, NX, nsym);
[ind_uAinuX, ind_notuAinuX] = getindwithin(NA, NX);

ind_AinX = [1; 1+(ind_usymAinusymX); 1+NsymX+(ind_uAinuX)]; 
ind_phiAinA = 1+NsymA+(1:NNA);


%% Precomputation of the tails (the \Delta^{-2} acting outside of Abar)
[~, usym, Phi] = extract_OK(X, NX, nsym, 'vectsym');

InvLap2phi_notAinX = InvLap2(ind_notuAinuX) .* Phi(ind_notuAinuX,:);

usym = usym(ind_notusymAinusymX,:);
InvLap2 = InvLap2(ind_usymXinuX);
InvLap2 = InvLap2(ind_notusymAinusymX);
InvLap2usym_notAinX = InvLap2 .* usym;

%% Y = A*X
% lambda part
Ylambda = Abar(1,:) * X(ind_AinX,:);

% usym part
Yusym = zeros(NsymX,size(X,2));
if useintervals
    Yusym = intval(Yusym);
end
Yusym(ind_usymAinusymX,:) = Abar(ind_usymAinA,:)*X(ind_AinX,:);
Yusym(ind_notusymAinusymX,:) = -InvLap2usym_notAinX;

% phi part
Yphi = zeros(NNX,size(X,2));
if useintervals
    Yphi = intval(Yphi);
end
Yphi(ind_uAinuX,:) = Abar(ind_phiAinA,:)*X(ind_AinX,:);
Yphi(ind_notuAinuX,:) = -InvLap2phi_notAinX;

Y = [Ylambda; 
     Yusym; 
     Yphi];
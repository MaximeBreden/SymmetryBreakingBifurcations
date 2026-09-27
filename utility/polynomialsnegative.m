function [rmin,eta] = polynomialsnegative(Y,Z,W)
% Tries to find the (approximately) minimal r = rmin 
% for which the radii polynomials 
% Y + (Z-I)*r + 1/2*r*W*r 
% as well as
% [Z-I+r*W]*eta 
% are negative
% Also returns eta in case you are interested
% If Y,Z,W are intervals the result is guaranteed via interval arithmetic
%
% Y is a column vector bounding the components ||T_m(x)-x_m|| <= Y(m)
% Z is a matrix bounding the derative components ||D_i T_m(x)|| <= Z(m,i)
% W is a tensor bounding ||D_i T_m(x+y) + D_i T_m(x)|| <=  sum_j W(m,i,j) ||y_j||
% If the W-bound holds for r<=r_* you will need to check that rmin<=r_* a posteriori
% i.e. this function does not check that.
% One may take W(i,j,k) a bound on sup_{|z|<=r_*} ||D_j D_k T_i(x+z)|| but not necessarily

M = length(Y);

% floating point versions
Yf = i2f(Y, 'sup');
Yf = max(Yf,realmin); %Yf == 0 causes unnecessary problems in the algorithm
Zf = i2f(Z, 'sup');
Wf = i2f(W, 'sup');
Wfreshape=reshape(Wf,[],M);
WWf =  @(r) reshape(Wfreshape*r,M,M); %sum over j assuming ||y|| <= r
Wfpermutereshape=reshape(permute(Wf,[1 3 2]),[],M); %transposed version needed for the derivative
WW2f =  @(r) reshape(Wfpermutereshape*r,M,M); 
Pf = @(r) Yf + Zf*r + WWf(r)*r/2 - r;
DPf = @(r) Zf + (WWf(r)+WW2f(r))/2 - eye(M);

% interval arithmetic versions
Wreshape=reshape(W,[],M);
WW =  @(r) reshape(Wreshape*r,M,M);
Pintval = @(r) Y + Z*r + WW(r)*r/2 -r;
% matrix
Matrix = @(r) Z + WW(r); 

% find approximate zero
r0 = zeros(M,1);
newtoncounter = 0;
maxnewtoniterates = 15;
newtontolerance = 1e-15;
convergencetolerance=1e-13;
dr0 = 1;
while newtoncounter <= maxnewtoniterates && max(abs(dr0)) > newtontolerance
    dr0 = -DPf(r0)\Pf(r0);
    if any(isnan(dr0))
        % there are some NaN 
        dr0(isnan(dr0))=0;
        nanflag=true;
    else
        nanflag=false;
    end
    r0 = r0+dr0;
    newtoncounter = newtoncounter+1;
end
if ~all(r0>0) || max(abs(dr0)) > convergencetolerance || nanflag
    disp('No good set of radii found for inclusion')
end
r0 = max(r0,realmin);

% determine search direction such that the components in Pf 
% decrease "equally" relative to the scale r0
direction = -DPf(r0)\r0;
% normalize the direction
direction = direction./r0;
direction = direction/max(direction);

% starting from the floating point zero of
% Y+(Z-I)*r+1/2*r*W*r 
% take a small step in direction to find a point
% where all these polynomials are negative
success = false;
partialsuccess = false;
rmin = NaN;
if all(r0 >= 0)
    r2 = r0; 
    n = -52;
    % search for negative point
    while n <= 150 && ~partialsuccess && all(r2>0)
        if all(Pintval(r2)<0) 
            rmin = r2;
            partialsuccess = true;
        else
            % step is search direction
            r2 = r0.*(ones(M,1)+direction*2^n);
            n = n+1;
        end
    end
end

if partialsuccess
    % check contractivity
    if M == 1
        if Matrix(rmin)-1 < 0
            success = true;
            eta = 1;
        end
    else
        [dominant,eta] = CollatzWielandt(Matrix(rmin));       
        if dominant-1 < 0
            success = true;
        else
            disp('Inclusion found, but no contraction.')
        end            
    end
end
         
if ~success 
    if M == 1
        disp('Radii polynomial not negative');
    else
        disp('Radii polynomials not negative simultaneously');
    end
    rmin = NaN;
    eta = NaN;
end

end

function [PFupperbound,testvector] = CollatzWielandt(A)
% provides upper bound on the spectral radius of A
% where A is assumed to have nonnegative elements
% (or irreducible non-negative matrix)
% The upper bound is rigorous in case A is an interval matrix.
%
% When A is non-negative but not necessarily irreducible
% (for example lower triangular)
% then the function returns an upper bound on 
% min_(x>0) max_i [Ax]i / xi 

% use some perturbations to make the algorithm more stable
perturbations=[0,1,10,100,1000]*eps;
upperbounds=zeros([length(perturbations),1]);
testvectors=zeros([size(A,2),length(perturbations)]);
useintervals = isa(A(1), 'intval');
if useintervals
    upperbounds = intval(upperbounds);
    testvectors = intval(testvectors);
end
for j=1:length(perturbations)
    Aeps=max(i2f(A,'mid'),perturbations(j));
    [eigenvectors,eigenvalues] = eig(Aeps);
    [~,m] = max(diag(eigenvalues));
    y = abs(eigenvectors(:,m));
    testvectors(:,j) = max(y,eps); % no division by 0
    % Collatz-Wielandt gives upper bound on Perron-Frobenius eigenvalue
    upperbounds(j) = max((A*testvectors(:,j))./testvectors(:,j));
end
% pick the best upperbound
[PFupperbound,jmin]=min(i2f(upperbounds,'sup'));
testvector=testvectors(:,jmin);

if useintervals
    PFupperbound = intval(PFupperbound);
end
end


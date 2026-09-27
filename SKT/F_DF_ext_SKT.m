function [Fext, DFext] = F_DF_ext_SKT(X, para, N, nsym, l)

% The extended system is called Fext here. This function computes Fext(X), 
% and its Frechet derivative DFext(X).

% X = (lambda; usym; phi)
[lambda, U1, U2, Phi1, Phi2] = extract_SKT(X, N, nsym, 'vectnosym');
U = [U1; U2];
Phi = [Phi1; Phi2];
para.d1 = lambda;
para.d2 = lambda;

if nargin < 5
    % If no linear map l is prescribed, we use l(v) = Phi^* v.
    l = Phi;
end

if nargout == 1
    [F, DuF] = F_DF_SKT(U, para, N);
else %if we also want DFext
    [F, DuF, DlambdaF, DlambdauFphi, DuuFphi] = F_DF_SKT(U, para, N, Phi);
end

[indsym, ~, Nsym] = getindsym(N, nsym);
NN = prod(N);
indsymU = [indsym; NN+indsym];

Fext = [sum(l.*Phi)-1;
        F(indsymU);
        DuF*Phi];

if nargout>1 %if we also want DFext
    DFext = [0, zeros(1,2*Nsym), l';
             DlambdaF(indsymU), DuF(indsymU,indsymU), zeros(2*Nsym,2*NN);
             DlambdauFphi, DuuFphi(:,indsymU), DuF];
end

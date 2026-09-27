function [Fext, DFext] = F_DF_ext_OK(X, para, N, nsym, l)

% The extended system is called Fext here. This function computes Fext(X), 
% and its Frechet derivative DFext(X).

% X = (lambda; usym; phi)
[lambda, U, Phi] = extract_OK(X, N, nsym, 'vectnosym');
para.lambda = lambda;

if nargin < 5
    l = Phi;
end

if nargout == 1
    [F, DuF] = F_DF_OK(U, para, N);
else %if we also want DFext
    [F, DuF, DlambdaF, DlambdauFphi, DuuFphi] = F_DF_OK(U, para, N, Phi);
end

[indsym, ~, Nsym] = getindsym(N, nsym);

Fext = [sum(l.*Phi)-1;
        F(indsym);
        DuF*Phi];

if nargout>1 %if we also want DFext
    NN = prod(N);
    DFext = [0, zeros(1,Nsym), l';
             DlambdaF(indsym), DuF(indsym,indsym), zeros(Nsym,NN);
             DlambdauFphi, DuuFphi(:,indsym), DuF];
end

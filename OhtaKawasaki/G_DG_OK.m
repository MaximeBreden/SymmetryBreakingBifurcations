function [G, DG] = G_DG_OK(Z, X, para, N, nsym)

% The map G corresponding to the eigenproblem, and its derivative.

[mu, Psi] = extract_OK_psi(Z, N, 'vect');
[lambda, U, Phi] = extract_OK(X, N, nsym, 'vectnosym');
para.lambda = lambda;

[~, DuF] = F_DF_OK(U, para, N);

G = [sum(Psi.*Phi)-1;
     DuF'*Psi - mu*Psi];
     
if nargout>1 %if we also want DG
    DG = [0, Phi';
          -Psi, DuF'-mu*eye(prod(N))];
end

function [G, DG] = G_DG_SKT(Z, X, para, N, nsym)

% The map G corresponding to the eigenproblem, and its derivative.

[mu, Psi1, Psi2] = extract_SKT_psi(Z, N, 'vect');
Psi = [Psi1; Psi2];

[lambda, U1, U2, Phi1, Phi2] = extract_SKT(X, N, nsym, 'vectnosym');
U = [U1; U2];
Phi = [Phi1; Phi2];

para.lambda = lambda;

[~, DuF] = F_DF_SKT(U, para, N);

G = [sum(Psi.*Phi)-1;
     DuF'*Psi - mu*Psi];
     
if nargout>1 %if we also want DG
    DG = [0, Phi';
          -Psi, DuF'-mu*eye(2*prod(N))];
end

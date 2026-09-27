function out = checkeig_OK(X, Z, rX, rZ,  N, nsym, nu, s)

% Checking whether psi^* phi \neq 0 (which proves that mu = 0). If we can
% prove that psi^* phi \neq 0, out = true, otherwise out = false.

% get Phi and Psi
[~, ~, Phi] = extract_OK(X, N, nsym);
[~, Psi] = extract_OK_psi(Z, N);

weights = getweights(N, nu, s-4);

normPhi = sum(abs(Phi).*weights);
normPsi = max(abs(Psi)./weights);

rPhi = rX(3);
rPsi = rZ(2);

out = ( normPsi*rPhi + normPhi*rPsi + rPhi*rPsi < abs(sum(Psi.*Phi)) );
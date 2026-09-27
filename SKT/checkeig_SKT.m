function out = checkeig_SKT(X, Z, rX, rZ,  N, nsym, nu, k)

% Checking whether psi^* phi \neq 0 (which proves that mu = 0). If we can
% prove that psi^* phi \neq 0, out = true, otherwise out = false.

% get Phi and Psi
[~, ~, ~, Phi1, Phi2] = extract_SKT(X, N, nsym);
[~, Psi1, Psi2] = extract_SKT_psi(Z, N);

weights = getweights(N, nu, k-2);

normPhi1 = sum(abs(Phi1).*weights);
normPhi2 = sum(abs(Phi2).*weights);
normPsi1 = max(abs(Psi1)./weights);
normPsi2 = max(abs(Psi2)./weights);

rPhi1 = rX(4);
rPhi2 = rX(5);
rPsi1 = rZ(2);
rPsi2 = rZ(3);

out = ( normPsi1*rPhi1 + normPsi2*rPhi2 + normPhi1*rPsi1 + normPhi2*rPsi2 + ...
        rPhi1*rPsi1 + rPhi2*rPsi2 < abs(sum(Psi1.*Phi1+Psi2.*Phi2)) );
function out = norm_blocks_SKT(X, nu, s, N, nsym)

% Takes the weighted \ell^1_k norm of each component of X

[lambda, U1, U2, Phi1, Phi2] = extract_SKT(X, N, nsym, 'vectnosym');

weights = getweights(N, nu, s)';

out = [abs(lambda);
       weights * abs(U1);
       weights * abs(U2);
       weights * abs(Phi1);
       weights * abs(Phi2)];
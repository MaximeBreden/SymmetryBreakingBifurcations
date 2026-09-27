function out = norm_blocks_OK(X, nu, s, N, nsym)

% Takes the weighted \ell^1_k norm of each component of X

[lambda, U, Phi] = extract_OK(X, N, nsym, 'vectnosym');

weights = getweights(N, nu, s)';

out = [abs(lambda);
       weights * abs(U);
       weights * abs(Phi)];
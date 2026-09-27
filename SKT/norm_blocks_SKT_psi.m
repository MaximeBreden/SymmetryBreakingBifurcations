function out = norm_blocks_SKT_psi(X, nu, s, N)

% Takes the weighted \ell^\infty_{-k} norm of each component of Z

[mu, Psi1, Psi2] = extract_SKT_psi(X, N, 'vect');

weights = getweights(N, nu, s);

out = [abs(mu);
       max(abs(Psi1)./weights);
       max(abs(Psi2)./weights)];
function out = norm_blocks_OK_psi(X, nu, s, N)

% Takes the weighted \ell^\infty_{-k} norm of each component of Z

[mu, Psi] = extract_OK_psi(X, N, 'vect');

out = [abs(mu);
       max(abs(Psi)./getweights(N, nu, s))];
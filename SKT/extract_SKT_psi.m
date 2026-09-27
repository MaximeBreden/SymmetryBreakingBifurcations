function [mu, Psi1, Psi2] = extract_SKT_psi(Z, N, type)

% Same thing as extract_SKT, but for Z = (mu; psi1; psi2) instead of 
% X = (lambda; u1; u2; phi1; phi2).

if nargin < 3
    type = 'vect';
end

NN = prod(N);
mu = Z(1,:);
Psi1 = Z(1+(1:NN),:);
Psi2 = Z(1+NN+(1:NN),:);

switch type
    case 'vect'
        % there is nothing more to do
    case 'mat'
        Psi1 = reshape(Psi1, N);
        Psi2 = reshape(Psi2, N);
end
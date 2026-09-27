function [mu, Psi] = extract_OK_psi(Z, N, type)

% Same thing as extract_OK, but for Z = (mu; psi) instead of 
% X = (lambda; u; phi).

if nargin < 3
    type = 'vect';
end

mu = Z(1,:);
Psi = Z(2:end,:);

switch type
    case 'vect'
        % there is nothing more to do
    case 'mat'
        Psi = reshape(Psi, N);
end
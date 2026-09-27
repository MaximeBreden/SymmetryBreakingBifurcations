function [lambda, U, Phi] = extract_OK(X, N, nsym, type)

% Extract the different components from X = (lambda; u; phi).
% In dimension 2, u and phi may be stored in matrices, or in vectors.
% Additionally, for u, one may want to only get the subset of coefficients 
% corresponding to the symmetry (which is what is stored in X), or all the 
% coefficients (with many zeros). These different options can be selected
% via the optionnal input 'type', which can be equal to 'mat', 'vectnosym'
% or 'vectsym' (default).

if nargin == 3
    type = 'vectsym';
end

[indsym, ~, Nsym] = getindsym(N, nsym);
NN = prod(N);

lambda = X(1,:);
U = X(1+(1:Nsym),:);
Phi = X(1+Nsym+(1:NN),:);

switch type
    case 'vectsym'
        % there is nothing more to do
    case 'vectnosym'
        u = zeros(N);
        if isa(X(1),'intval')
            u = intval(u);
        end
        u(indsym) = U;
        U = u(:);
    case 'mat'
        u = zeros(N);
        if isa(X(1),'intval')
            u = intval(u);
        end
        u(indsym) = U;
        U = u;
        Phi = reshape(Phi, N);
end
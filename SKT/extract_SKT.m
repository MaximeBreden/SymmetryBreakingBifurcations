function [lambda, U1, U2, Phi1, Phi2] = extract_SKT(X, N, nsym, type)

% Extract the different components from X = (lambda; u1; u2; phi1; phi2).
% In dimension 2, ui and phii may be stored in matrices, or in vectors.
% Additionally, for ui, one may want to only get the subset of coefficients 
% corresponding to the symmetry (which is what is stored in X), or all the 
% coefficients (with many zeros). These different options can be selected
% via the optionnal input 'type', which can be equal to 'mat', 'vectnosym'
% or 'vectsym' (default).

if nargin == 3
    type = 'vectsym';
end

[indsym, ~, Nsym] = getindsym(N, nsym); % indices corresponding the the symmetry
NN = prod(N);

lambda = X(1,:);
U1 = X(1+(1:Nsym),:);
U2 = X(1+Nsym+(1:Nsym),:);
Phi1 = X(1+2*Nsym+(1:NN),:);
Phi2 = X(1+2*Nsym+NN+(1:NN),:);

switch type
    case 'vectsym'
        % there is nothing more to do
    case 'vectnosym'
        u1 = zeros(N);
        u2 = zeros(N);
        if isa(X(1),'intval')
            u1 = intval(u1);
            u2 = intval(u2);
        end
        u1(indsym) = U1;
        u2(indsym) = U2;
        U1 = u1(:);
        U2 = u2(:);
    case 'mat'
        u1 = zeros(N);
        u2 = zeros(N);
        if isa(X(1),'intval')
            u1 = intval(u1);
            u2 = intval(u2);
        end
        u1(indsym) = U1;
        u2(indsym) = U2;
        U1 = u1;
        U2 = u2;
        Phi1 = reshape(Phi1, N);
        Phi2 = reshape(Phi2, N);
end
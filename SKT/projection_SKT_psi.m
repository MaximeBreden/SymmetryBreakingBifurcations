function Znew = projection_SKT_psi(Z, Nnew, N)

% Change the size psi to Nproj (pad by zeros or truncates)

[mu, psi1, psi2] = extract_SKT_psi(Z, N, 'mat');

% psi
NN = min(N, Nnew);
psi1new = zeros(Nnew);
psi2new = zeros(Nnew);
if isa(Z(1),'intval')
    psi1new = intval(psi1new);
    psi2new = intval(psi2new);
end
psi1new(1:NN(1),1:NN(2)) = psi1(1:NN(1),1:NN(2));
psi2new(1:NN(1),1:NN(2)) = psi2(1:NN(1),1:NN(2));

% Z
Znew = [mu; psi1new(:); psi2new(:)];



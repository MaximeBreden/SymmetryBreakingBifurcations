function Znew = projection_OK_psi(Z, Nnew, N)

% Change the size psi to Nproj (pad by zeros or truncates)

[mu, psi] = extract_OK_psi(Z, N, 'mat');

% psi
NN = min(N, Nnew);
psinew = zeros(Nnew);
if isa(Z(1),'intval')
    psinew = intval(psinew);
end
psinew(1:NN(1),1:NN(2)) = psi(1:NN(1),1:NN(2));

% Z
Znew = [mu; psinew(:)];



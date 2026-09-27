function Xnew = projection_OK(X, Nnew, N, nsym)

% Change the size of usym and phi to Nproj (pad by zeros or truncates)

[lambda, u, phi] = extract_OK(X, N, nsym, 'mat');

% phi
NN = min(N, Nnew);
phinew = zeros(Nnew);
if isa(X(1),'intval')
    phinew = intval(phinew);
end
phinew(1:NN(1),1:NN(2)) = phi(1:NN(1),1:NN(2));

% u
unew = zeros(Nnew);
if isa(X(1),'intval')
    unew = intval(unew);
end
unew(1:NN(1),1:NN(2)) = u(1:NN(1),1:NN(2));

indsymnew = getindsym(Nnew, nsym);
Usymnew = unew(indsymnew);

% X
Xnew = [lambda; Usymnew; phinew(:)];



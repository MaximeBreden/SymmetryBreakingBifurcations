function Xnew = projection_SKT(X, Nnew, N, nsym)

% Change the size of usym and phi to Nproj (pad by zeros or truncates)

[lambda, u1, u2, phi1, phi2] = extract_SKT(X, N, nsym, 'mat');

% phi
NN = min(N, Nnew);
phi1new = zeros(Nnew);
phi2new = zeros(Nnew);
if isa(X(1),'intval')
    phi1new = intval(phi1new);
    phi2new = intval(phi2new);
end
phi1new(1:NN(1),1:NN(2)) = phi1(1:NN(1),1:NN(2));
phi2new(1:NN(1),1:NN(2)) = phi2(1:NN(1),1:NN(2));

% u
u1new = zeros(Nnew);
u2new = zeros(Nnew);
if isa(X(1),'intval')
    u1new = intval(u1new);
    u2new = intval(u2new);
end
u1new(1:NN(1),1:NN(2)) = u1(1:NN(1),1:NN(2));
u2new(1:NN(1),1:NN(2)) = u2(1:NN(1),1:NN(2));

indsymnew = getindsym(Nnew, nsym);
U1symnew = u1new(indsymnew);
U2symnew = u2new(indsymnew);

% X
Xnew = [lambda; U1symnew; U2symnew; phi1new(:); phi2new(:)];



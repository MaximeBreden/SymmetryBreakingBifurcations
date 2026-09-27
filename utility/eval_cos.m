function ueval = eval_cos(u, x, y)

% u is a matrix of cosine coefficients

[Kx, Ky] = size(u);
kx = 0 : Kx-1;
ky = 0 : Ky-1;

Mx = cos(x(:)*kx);
My = cos(y(:)*ky);

u(2:end,:) = sqrt(2)*u(2:end,:);
u(:,2:end) = sqrt(2)*u(:,2:end);

ueval = reshape( sum ( (Mx*u).*My, 2), size(x));
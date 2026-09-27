function w = getweights(N, nu, s)

% Generates the weights \xi_n^{(s)} * \nu^|n|, with \xi_n^{(s)} as in the
% paper.

if isa(nu, 'intval')
    sqrt2 = sqrt(intval(2));
    N1 = intval((0:N(1)-1)');
    N2 = intval((0:N(2)-1));
else
    sqrt2 = sqrt(2);
    N1 = (0:N(1)-1)';
    N2 = (0:N(2)-1);
end

% the part coming from the scalar product (sqrt(2) for n>=1, in each direction)
w = 2*ones(N);
if isa(nu, 'intval')
    w = intval(w);
end
w(1,1) = 1;
w(1,2:end) = sqrt2;
w(2:end,1) = sqrt2;

% the geometric part (nu is always taken equal to 1 in the paper)
w = w .* nu.^(N1+N2);

% the algebraic part
w = w .* (1+sqrt(N1.^2+N2.^2)).^s;

% reshaping into a vector
w = w(:);
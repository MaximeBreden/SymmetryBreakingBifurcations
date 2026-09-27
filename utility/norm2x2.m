function out = norm2x2(w, nu, s)

out = zeros(2,2);
if isa(nu,'intval')
    out = intval(out);
end

for i = 1:2
    for j = 1:2
        u = w{i,j};
        N = size(u);
        out(i,j) = sum(getweights(N, nu, s) .* abs(u(:)), "all");
    end
end
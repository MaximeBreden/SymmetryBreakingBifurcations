function out = opnorm_blocks_OK_psi(B, nu, sr, sc, Nr, Nc)

% Takes the weighted \ell^1 operator norm of each block, from \ell^1_{sc}
% to \ell^1_{sr}

weightsr = getweights(Nr, nu, sr);
weightsc = getweights(Nc, nu, sc);

NNr = prod(Nr);
NNc = prod(Nc);

indr2 = 1+(1:NNr);
indc2 = 1+(1:NNc);

B = abs(B);

out = zeros(2,2);
if isa(B(1), 'intval')
    out = intval(out);
end

% Row 1
out(1,1) = B(1,1);
out(1,2) = max(B(1,indc2)./weightsc');

% Row 2
out(2,1) = weightsr'*B(indr2,1);
out(2,2) = max((weightsr'*B(indr2,indc2))./weightsc');



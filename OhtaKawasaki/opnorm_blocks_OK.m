function out = opnorm_blocks_OK(B, nu, sr, sc, Nr, Nc, nsym)

% Takes the weighted \ell^1 operator norm of each block, from \ell^1_{sc}
% to \ell^1_{sr}

[indsymr, ~, Nsymr] = getindsym(Nr, nsym);
[indsymc, ~, Nsymc] = getindsym(Nc, nsym);

weightsr = getweights(Nr, nu, sr)';
weightsc = getweights(Nc, nu, sc)';

weightssymr = weightsr(indsymr);
weightssymc = weightsc(indsymc);

NNr = prod(Nr);
NNc = prod(Nc);

indr2 = 1+(1:Nsymr);
indr3 = 1+Nsymr+(1:NNr);
indc2 = 1+(1:Nsymc);
indc3 = 1+Nsymc+(1:NNc);

B = abs(B);

out = zeros(3,3);
if isa(B(1), 'intval')
    out = intval(out);
end

% Row 1
out(1,1) = B(1,1);
out(1,2) = max(B(1,indc2)./weightssymc);
out(1,3) = max(B(1,indc3)./weightsc);

% Row 2
out(2,1) = weightssymr*B(indr2,1);
out(2,2) = max((weightssymr*B(indr2,indc2))./weightssymc);
out(2,3) = max((weightssymr*B(indr2,indc3))./weightsc);

% Row 3
out(3,1) = weightsr*B(indr3,1);
out(3,2) = max((weightsr*B(indr3,indc2))./weightssymc);
out(3,3) = max((weightsr*B(indr3,indc3))./weightsc);


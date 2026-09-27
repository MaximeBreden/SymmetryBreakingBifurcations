function [iX, ipara, iAbar] = converttointval_OK(X, para, Abar)

% Turns everything into intervals. 
% /!\ For the examples from the paper this is fine, because sigma and Lx
% are always representable exactly as floating point numbers, but otherwise
% one has to be slightly more careful.

iX = intval(X);

ipara.lambda = intval(para.lambda);
ipara.mu = midrad(para.mu, eps(para.mu));
ipara.sigma = intval(para.sigma);
ipara.Lx = intval(para.Lx);
if isfield(para,'Ly')
    ipara.Ly = midrad(para.Ly, eps(para.Ly));
end

iAbar = intval(Abar);

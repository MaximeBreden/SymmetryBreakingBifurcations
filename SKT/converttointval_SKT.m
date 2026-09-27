function [iX, ipara, iAbar, iw, igamma] = converttointval_SKT(X, para, Abar, w, gamma)

% Turns everything into intervals. 
% /!\ For the examples from the paper this is fine, because all the
% parameters of the model are representable exactly as floating point
% numbers, but otherwise one has to be slightly more careful.

iX = intval(X);

ipara.d1 = iX(1);
ipara.d2 = iX(1);
ipara.d11 = intval(para.d11);
ipara.d12 = intval(para.d12);
ipara.d21 = intval(para.d21);
ipara.d22 = intval(para.d22);
ipara.r1 = intval(para.r1);
ipara.a1 = intval(para.a1);
ipara.b1 = intval(para.b1);
ipara.r2 = intval(para.r2);
ipara.a2 = intval(para.a2);
ipara.b2 = intval(para.b2);
ipara.Lx = intval(para.Lx);
if isfield(para,'Ly')
    ipara.Ly = intval(para.Ly);
end

iAbar = intval(Abar);

iw = {intval(w{1,1}), intval(w{1,2});
      intval(w{2,1}), intval(w{2,2})};

igamma = {intval(gamma{1,1}), intval(gamma{1,2});
          intval(gamma{2,1}), intval(gamma{2,2})};
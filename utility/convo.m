function w = convo(u, v, Nw)
%Compute the convolution product of two vectors of cosine coeffs

%The two vector can be of different lengths.

%The third input prescribes the length of the output and is optional. 

Nv = size(v);

% V = v(:);
% W = convomat(u, Nw, Nv) * V;
% w = reshape(W, Nw);

w = reshape(convomat(u, Nw, Nv) * v(:), Nw);
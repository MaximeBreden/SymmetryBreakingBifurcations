function [indsym, indcompl, Nsym] = getindsym(N, nsym)

% Computes the indices indsym so that uisym = ui(indsym). The remaining
% indices are given in indcompl, and Nsym gives the numbers of elements
% remaining in uisym.

% EXAMPLE. With nsym = [1,2] and N = [3,6], the symmetric coefficients are 
% the ones indicates by a 1 in the following array
% [1 0 0 0 1 0
%  0 0 1 0 0 0
%  1 0 0 0 1 0],
% i.e., they correspond to the indices [1,3,8,13,15] (with linear indexing).
% Thus, indsym = [1,3,8,13,15], 
% incdompl = [2,4,5,6,7,9,10,11,12,14,16,17,18], and Nsym = 5.

ind = reshape(1:(N(1)*N(2)), N);

mask = false(N);
if isscalar(nsym) % 1D case
    mask(1:nsym:N(1)) = true;
else % 2D case   
    mask(1:2*nsym(1):N(1),1:2*nsym(2):N(2)) = true;
    mask(nsym(1)+1:2*nsym(1):N(1),nsym(2)+1:2*nsym(2):N(2)) = true;   
end

indsym = ind(mask);
indcompl = ind(not(mask));
Nsym = length(indsym);
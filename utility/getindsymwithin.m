function [indsymwithin, indsymoutside, indoutsidesym, Nsym] = getindsymwithin(Nsmall, N, nsym)

% Computes a bunch of index sets, related to having a set of "symmetric 
% coefficients" in an element of size N, and wanting to keep those for a
% smaller size. This is better explained on an example.

% EXAMPLE. With nsym = [1,2] and N = [3,6], the symmetric coefficients are 
% the ones indicates by 1 in the following array
% [1 0 0 0 1 0
%  0 0 1 0 0 0
%  1 0 0 0 1 0],
% i.e., they correspond to the indices [1,3,8,13,15] (with linear indexing).
% If we restrict ourselves to Nsmall = [2,3], we want to keep only the 1's 
% indicated below
% [1 0 0 0 0 0
%  0 0 1 0 0 0
%  0 0 0 0 0 0],
% i.e., the indices [1,8].
%
% The first output indsymwithin gives the position of the "symmetric 
% coefficients" for Nsmall within those for N. In this example 
% indsymwithin = [1,3] (because 1 and 8 are in position 1 and 3 in
% [1,3,8,13,15]).
%
% The second output indsymoutside gives the position of the remaining
% "symmetric coefficients" within those for N. In this example 
% indsymoutside = [2,4,5] (corresponding to 3, 13 and 15).
%
% The third output indoutsidesym also gives the position of the remaining
% "symmetric coefficients", but within all the coefficients (and not only
% within the symetric ones). For this example, indoutsidesym = [3, 13, 15].
%
% Nsym is simply the total number of "symmetric coefficients", so 5 here.

mask = false(N);
maskwithin = false(N);
if isscalar(nsym) % 1D case
    mask(1:nsym:N(1)) = true;
    maskwithin(1:nsym:Nsmall(1)) = true;
    maskoutside = mask;
    maskoutside(1:nsym:Nsmall(1)) = false;
else % 2D case   
    mask(1:2*nsym(1):N(1),1:2*nsym(2):N(2)) = true;
    mask(nsym(1)+1:2*nsym(1):N(1),nsym(2)+1:2*nsym(2):N(2)) = true;

    maskwithin(1:2*nsym(1):Nsmall(1),1:2*nsym(2):Nsmall(2)) = true;
    maskwithin(nsym(1)+1:2*nsym(1):Nsmall(1),nsym(2)+1:2*nsym(2):Nsmall(2)) = true;

    maskoutside = mask;
    maskoutside(1:2*nsym(1):Nsmall(1),1:2*nsym(2):Nsmall(2)) = false;
    maskoutside(nsym(1)+1:2*nsym(1):Nsmall(1),nsym(2)+1:2*nsym(2):Nsmall(2)) = false;
end

ind = reshape(1:(N(1)*N(2)), N);
indsym = ind(mask);
Nsym = length(indsym);

indlin = zeros(N);
indlin(mask) = 1:Nsym;

indsymwithin = indlin(maskwithin);
indsymoutside = indlin(maskoutside);
indoutsidesym = ind(maskoutside);
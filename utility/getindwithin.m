function [indwithin, indoutside] = getindwithin(Nsmall, N)

% Computes the indices describing an element of size Nsmall within a
% larger element of size N

% EXAMPLE. With Nsmall = [2,3] and N = [3,6], we want to extract the entries 
% indicates by a 1 in the following array
% [1 1 1 0 0 0
%  1 1 1 0 0 0
%  0 0 0 0 0 0].
% Thus, indwithin = [1,2,4,5,7,8], and the remaining indices are
% indoutside = [3,6,9,10,11,12,13,14,15,16,17,18].

ind = reshape(1:(N(1)*N(2)), N);
maskwithin = false(N);
maskwithin(1:Nsmall(1),1:Nsmall(2)) = true;
indwithin = ind(maskwithin);
indoutside = ind(not(maskwithin));

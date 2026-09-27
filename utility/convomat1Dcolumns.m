function M = convomat1Dcolumns(u, Noutput, Ninput)
%Compute the matrix associated to the convolution product of a vector of 
%cosine coeffs with a vector of cosine coeffs. That is, u is assumed to be
%a vector of cosine coeffs, and M is such that, for any vector v of cosine
%coeffs, u*v=Mv.
%
%The assumed normalizations is u=u_0+2*\sum_{k\geq 1} u_k cos(kx) 
%
%The second and third input are optional, and can be used to enforce the 
%size of M. By default, the size is the one of u.

N2 = size(u,2);
M = zeros(Noutput, Ninput*N2);
if isa(u(1),'intval')
    M = intval(M);
end
for n = 1:N2
    M(:,(n-1)*Ninput+(1:Ninput)) = convomat1D(u(:,n), Noutput, Ninput);
end


function M = convomat1D(u, Noutput, Ninput)
%Compute the matrix associated to the convolution product of a vector of 
%cosine coeffs with a vector of cosine coeffs. That is, u is assumed to be
%a vector of cosine coeffs, and M is such that, for any vector v of cosine
%coeffs, u*v=Mv.
%
%The assumed normalizations is u=u_0+2*\sum_{k\geq 1} u_k cos(kx) 
%
%The second and third input are optional, and can be used to enforce the 
%size of M. By default, the size is the one of u.

%% Padding u by zeros if needed for the output size
if nargin >= 2
    N = length(u);
    if nargin == 2
        Ninput = Noutput;
    end
    Next = max(Noutput, Ninput);
    if exist('intval','file') && isintval(u(1))
        u = [u;intval(zeros(Next-N,1))];
    else
        u = [u;zeros(Next-N,1)];
    end
end

%% Construction of the multiplication operator M
A = hankel(u);
C = toeplitz(u,u); % /!\ If u is complex, toeplitz(u) gives toeplitz(u,u'), which is not what we want here !!
C(:,1) = 0;
M = A+C;

%% Truncating M
if nargin >= 2 
    M = M(1:Noutput,1:Ninput);
end

function M = convomat(u, Noutput, Ninput)
%Compute the matrix associated to the convolution product of a vector of 
%cosine coeffs with a vector of cosine coeffs. That is, u is assumed to be
%a vector of cosine coeffs, and M is such that, for any vector v of cosine
%coeffs, u*v=Mv.
%
%The assumed normalizations is u=u_0+sqrt(2)*\sum_{k\geq 1} u_k cos(kx) 
%M acts on elements of size Nintput = N, and then yield and element of size Noutput = N
%The second and third input are optional, and can be used to enforce the 
%size of M. By default, Ninput = N and Noutput = N

%% Default sizes
N = size(u);
if nargin < 3
    Ninput = N;
    if nargin < 2
        Noutput = N;
    end
end

%% Rescaling of u
if isa(u(1), 'intval')
    sqrt2 = sqrt(intval(2));
else
    sqrt2 = sqrt(2);
end

resc_u = 2*ones(N);
if isa(u(1), 'intval')
    resc_u = intval(resc_u);
end
resc_u(1,1) = 1;
resc_u(1,2:end) = sqrt2;
resc_u(2:end,1) = sqrt2;
u = u./resc_u;

%% Construction of the multiplication operator M
H = hankel((1:Noutput(2))', [Noutput(2) zeros(1,Ninput(2)-1)]);
T = toeplitz((1:Noutput(2))', (1:Ninput(2))');
T(:,1) = 0;

% Largest index needed in the second direction
N2 = max(Noutput(2),Ninput(2));

% Adding row(s) of zeros in u (for the zeros in the Hankel matrix)
u = [zeros(N(1),1), u, zeros(N(1),N2-N(2))];

% Generating all the 1D convolution matrices (in a single matrix)
M = convomat1Dcolumns(u, Noutput(1), Ninput(1));

if not(isa(u(1),'intval'))
    % Extracting all the 1D convolution matrices separately 
    M = mat2cell(M, Noutput(1), Ninput(1)*ones(1,N2+1));
    
    % Generating the 2D convolution matrix
    M = cell2mat(M(H+1)) + cell2mat(M(T+1));
else 
    % We do the same, but going via midpoint and radius to speed up the
    % code (directly applying mat2cell and cell2mat to intvals works but is
    % somehow much slower).

    % Extracting all the 1D convolution matrices separately (also
    % separating midpoint and radius)
    Mmid = mat2cell(mid(M), Noutput(1), Ninput(1)*ones(1,N2+1));
    Mrad = mat2cell(rad(M), Noutput(1), Ninput(1)*ones(1,N2+1));

    % Generating the 2D convolution matrix
    MmidH = cell2mat(Mmid(H+1));
    MradH = cell2mat(Mrad(H+1));
    MmidT = cell2mat(Mmid(T+1));
    MradT = cell2mat(Mrad(T+1));
    M = midrad(MmidH,MradH) + midrad(MmidT,MradT);
end


%% Rescaling of M
resc_input = 2*ones(Ninput);
if isa(u(1), 'intval')
    resc_input = intval(resc_input);
end
resc_input(1,1) = 1;
resc_input(1,2:end) = sqrt2;
resc_input(2:end,1) = sqrt2;

resc_output = 2*ones(Noutput);
if isa(u(1), 'intval')
    resc_output = intval(resc_output);
end
resc_output(1,1) = 1;
resc_output(1,2:end) = sqrt2;
resc_output(2:end,1) = sqrt2;

M = resc_output(:) .* M ./ resc_input(:)';
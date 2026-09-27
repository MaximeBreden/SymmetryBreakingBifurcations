%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FILE TO RUN THE SKT PROOFS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% Change filename below to select the example you want from the paper.
% We have SKT1D_i for i = 1 to 4 (lines 1 to 4 of Table 2), 
% and SKT2D_i for i = 1,2 (lines 5 and 6 of Table 2).
filename = 'SKT1D_1'; 
load(['data/',filename,'.mat'], 'X', 'Z', 'para', 'N', 'nsym')

%% Choice of space (i.e., of weight for the norm)
% The weights behave like n^s. Any s>=0 is legit for proving the existence of a bifurcation.
% If you have a transcritical bifurcation and want to rigorously prove it, you need to take s>=2.
s = 0; % /!\ Take s = 2 if you care about proving that the bifurcation is transcritical

% Extra geometric weight nu^n (nu>=1). Since we ended up always taking
% nu=1, this was not even introduced in the paper.
nu = 1; 
if exist('intval','file')
    inu = intval('1');
end

%% Plot
[lambda, u1, u2, phi1, phi2] = extract_SKT(X, N, nsym, 'mat');
plotdata_SKT(u1, u2, phi1, phi2, para)
d = lambda;
v0 = sum(u2(:).*getweights(N,1,0));
normu = sqrt(sum(u1(:).^2));
sgtitle(['$d$ = ',num2str(d),', $v(0)$ = ',num2str(v0),', $\Vert u\Vert_{L^2}$ = ',num2str(normu)], 'FontSize', 14, 'Interpreter', 'latex');
drawnow
% saveas(gcf,[filename,'.eps'],'epsc')

%% "Prevalidation" (without interval arithmetic)
[rX, Abar, w, gamma] = proof_SKT(X, para, N, nsym, nu, s);
succesX = not(any(isnan(rX)));
if succesX
    fprintf("\nPrevalidation successful, with error bounds:\n");
    disp(rX)
end

%%  Rigorous proof (with interval arithmetic), you need Intlab for this
if exist('intval','file') && succesX
    [iX, ipara, iAbar, iw, igamma] = converttointval_SKT(X, para, Abar, w, gamma);
    irX = proof_SKT(iX, ipara, N, nsym, inu, s, iAbar, iw, igamma);
    isuccesX = not(any(isnan(irX)));
    if isuccesX
        fprintf("\nValidation successful, with error bounds:\n");
        disp(irX)
    end
end

%% This should be zero iff the bifurcation is a pitchfork
testbif = bifurcationtest_SKT(X, Z, para, N, nsym);
fprintf("\nApproximation of the quantity that should be zero iff the bifurcation is a pitchfork: %e", testbif)
if abs(testbif) < 1e-5
    fprintf("\nIt looks like the bifurcation might be a pitchfork (you need to prove this by hand).\n")
    return
else
    fprintf("\nIt looks like the bifurcation might be transcritical. ")
    if succesX
        fprintf("We now try to prove it.\n")
    else
        return
    end
end

%% Enclosure of Psi and rigorous determination of the bifurcation type

% "Prevalidation" (without interval arithmetic)
transcritical = false; % will become true only if we can prove that the bifurcation is indeed transcritical
[rZ, AbarG] = proof_SKT_psi(Z, X, para, N, nsym, nu, s, rX, w);
if any(isnan(rZ))
    fprintf("\nUnable to validate the approximate left eigenpair of DuF\n")
else
    fprintf("\nSuccesful prevalidation of a left eigenpair of DuF\n")
    if checkeig_SKT(X, Z, rX, rZ, N, nsym, nu, s)
        fprintf("\nWe have validated the correct Psi\n")
        [testbif, err] = bifurcationtest_SKT(X, Z, para, N, nsym, nu, s, rX, rZ);
        if err < abs(testbif)
            transcritical = true;
            fprintf("\nWe have proven, up to interval arithmetic, that the bifurcation is indeed transcritical\n")
        else
            fprintf("\nUnable to rigorously determine the type of the bifurcation\n")
        end
    else
        fprintf("\nUnable to prove that mu is 0\n")
    end
end

%  Rigorous proof (with interval arithmetic), you need Intlab for this
if exist('intval','file') && transcritical && isuccesX
    iZ = intval(Z);
    iAbarG = intval(AbarG);
    irZ = proof_SKT_psi(iZ, iX, ipara, N, nsym, inu, s, irX, iw, iAbarG);
    if any(isnan(irZ))
        fprintf("\nUnable to validate the approximate left eigenpair of DuF\n")
    else
        fprintf("\nSuccesful rigorous validation of a left eigenpair of DuF\n")
        if checkeig_SKT(iX, iZ, irX, irZ, N, nsym, inu, s)
            fprintf("\nWe have validated the correct Psi\n")
            [itestbif, ierr] = bifurcationtest_SKT(iX, iZ, ipara, N, nsym, inu, s, irX, irZ);
            if ierr < abs(itestbif)
                fprintf("\nWe have proven that the bifurcation is indeed transcritical\n")
            else
                fprintf("\nUnable to rigorously determine the type of the bifurcation\n")
            end
        else
            fprintf("\nUnable to prove that mu is 0\n")
        end
    end
end
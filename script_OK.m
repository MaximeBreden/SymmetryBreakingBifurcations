%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FILE TO RUN THE OHTA-KAWASAKI PROOFS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% Change filename below to select the example you want from the paper.
% We have OK1D_i for i = 1,2 (lines 1 and 2 of Table 1), 
% and OK2D_i for i = 1,2 (lines 3 and 4 of Table 1).
filename = 'OK1D_1'; 
load(['data/',filename,'.mat'], 'X', 'Z', 'para', 'N', 'nsym')

%% Choice of space (i.e., of weight for the norm)
% The weights behave like n^s. Any s>=0 is legit for proving the existence of a bifurcation.
% If you have a transcritical bifurcation and want to rigorously prove it, you need to take s>=2.
s = 0; % /!\ Take s = 2 if you care about proving that the bifurcation is transcritical

% Extra geometric weight nu^n (nu>=1). Since we ended up always taking
% nu=1, this was not even introduced in the paper.
nu = 1; 

% a priori radius for the proof
rstar = 1e-4;

if exist('intval','file')
    inu = intval('1');
    irstar = intval(rstar);
end

%% Plot
[lambda, u, phi] = extract_OK(X, N, nsym, 'mat');
plotdata_OK(u, phi, para)
u0 = sum(u(:).*getweights(N,1,0));
normu = sqrt(sum(u(:).^2));
sgtitle(['$\lambda$ = ',num2str(lambda),', $u(0)$ = ',num2str(u0),', $\Vert u\Vert_{L^2}$ = ',num2str(normu)], 'FontSize', 14, 'Interpreter', 'latex');
drawnow
% saveas(gcf,[filename,'.eps'],'epsc')


%% "Prevalidation" (without interval arithmetic)
[rX, Abar] = proof_OK(X, para, N, nsym, nu, s, rstar);
succesX = not(any(isnan(rX)));
if succesX
    fprintf("\nPrevalidation successful, with error bounds:\n");
    disp(rX)
end

%%  Rigorous proof (with interval arithmetic), you need Intlab for this
if exist('intval','file') && succesX
    [iX, ipara, iAbar] = converttointval_OK(X, para, Abar);
    irX = proof_OK(iX, ipara, N, nsym, inu, s, irstar, iAbar);
    isuccesX = not(any(isnan(irX)));
    if isuccesX
        fprintf("\nValidation successful, with error bounds:\n");
        disp(irX)
    end
end

%% This should be zero iff the bifurcation is a pitchfork
testbif = bifurcationtest_OK(X, Z, para, N, nsym);
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
[rZ, AbarG] = proof_OK_psi(Z, X, para, N, nsym, nu, s, rX);
if any(isnan(rZ))
    fprintf("\nUnable to validate the approximate left eigenpair of DuF\n")
else
    fprintf("\nSuccesful prevalidation of a left eigenpair of DuF\n")
    if checkeig_OK(X, Z, rX, rZ, N, nsym, nu, s)
        fprintf("\nWe have validated the correct Psi\n")
        [testbif, err] = bifurcationtest_OK(X, Z, para, N, nsym, nu, s, rX, rZ);
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
    irZ = proof_OK_psi(iZ, iX, ipara, N, nsym, inu, s, irX, iAbarG);
    if any(isnan(irZ))
        fprintf("\nUnable to validate the approximate left eigenpair of DuF\n")
    else
        fprintf("\nSuccesful rigorous validation of a left eigenpair of DuF\n")
        if checkeig_OK(iX, iZ, irX, irZ, N, nsym, inu, s)
            fprintf("\nWe have validated the correct Psi\n")
            [itestbif, ierr] = bifurcationtest_OK(iX, iZ, ipara, N, nsym, inu, s, irX, irZ);
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
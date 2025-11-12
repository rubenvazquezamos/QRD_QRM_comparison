% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.

% This code runs a numerical simulation of an acoustic 
% diffuser using COMSOL® and compares the obtained
% diffusion coefficient to that obtained using the
% transfer matrix method (TMM).

% Based on SchroederDiffuser_COMSOL

run("config.m)");

%       FEM MODELLING
%-------------------------------------------------------------------------%
save_dlg = false;
% Set save_dlg to true if you want to have the option to clear mesh and solutions 
% data after running the model
tStart = tic;
Ps_1 = QR_5(Freq,Geo,Probe,File,save_dlg); %call COMSOL model for QRD
Ps_2 = QRM(Freq,Geo,Probe,File,save_dlg);
Psflatnum = flat_plane(Freq,Geo,Probe,File,save_dlg);  %call COMSOL model for flat plane
tEnd = toc(tStart);
fprintf('FEM. time: %d minutes and  %.f seconds\n', floor(tEnd/60), rem(tEnd,60));
%-------------------------------------------------------------------------%

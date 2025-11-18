%% COMSOL PROBE INFORMATION
%-------------------------------------------------------------------------%
Probe.domain = 4; %radius of air domain in meters
Probe.radius = 3; %radius of arc in meters
Probe.theta_min =  0; %arc starting angle
Probe.theta_max = pi; %arc end angle
Probe.Resolution = 181;
Probe.theta_vector = linspace(Probe.theta_min,Probe.theta_max,Probe.Resolution);
Probe.Coordinates(1,:) = Probe.radius*cos(Probe.theta_vector); %x coordinates
Probe.Coordinates(2,:) = Probe.radius*sin(Probe.theta_vector); %y coordinates

%% COMSOL FILE INFORMATION
%-------------------------------------------------------------------------%
File.Path = [pwd,filesep,'Models'];
File.Tag1 = 'Comsol_QRD5';
File.Tag2 = 'flat_plane';
File.Tag3 = 'optgeo.mphbin';
File.Tag4 = 'QRM_5';
File.Extension = '.mph';

%       FEM MODELLING
%-------------------------------------------------------------------------%
save_dlg = false;
% Set save_dlg to true if you want to have the option to clear mesh and solutions 
% data after running the model
tStart = tic;
Ps.QRD_COMSOL = QR_5(Freq,Geo,Probe,File,save_dlg);
Ps.QRM_COMSOL = QRM(Freq,Geo,Probe,File,save_dlg);
Ps.flat_COMSOL = flat_plane(Freq,Geo,Probe,File,save_dlg);
tEnd = toc(tStart);
fprintf('FEM. time: %d minutes and  %.f seconds\n', floor(tEnd/60), rem(tEnd,60));
%-------------------------------------------------------------------------%
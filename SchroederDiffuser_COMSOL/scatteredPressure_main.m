% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.

% This code runs a numerical simulation of an acoustic 
% diffuser using COMSOL® and compares the obtained
% diffusion coefficient to that obtained using the
% transfer matrix method (TMM).

% Based on SchroederDiffuser_COMSOL

%% SESSION START UP COMMANDS
%-------------------------------------------------------------------------%
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',18)
set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultFigureColor',[1 1 1])
path = convertCharsToStrings(fileparts(matlab.desktop.editor.getActiveFilename));
cd(path)

clear; clear global *; clc; warning off; close all;

%% COMSOL FILE INFORMATION
%-------------------------------------------------------------------------%
File.Path = [pwd,filesep,'Models'];
File.Tag1 = 'Comsol_QRD5';
File.Tag2 = 'flat_plane';
File.Tag3 = 'optgeo.mphbin';
File.Extension = '.mph';
%-------------------------------------------------------------------------%

%% FREQUENCY
%-------------------------------------------------------------------------%
Freq.f_min = 250;                                  % Minimum Freq of interest
Freq.f_max = 4000;                               % Maximum Freq of interest
Freq.df = 50;                                    % Freq discretization
Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);    % Freq vector
Freq.Nf = numel(Freq.Vector);                   % Number of frequencies

%% GEOMETRY
%-------------------------------------------------------------------------%
N = 5; %prime number used to generate s_n
Geo.W = 7e-2; %width of single well
Geo.Li = 1e-2; %width in between wells
Geo.D = N*(Geo.W+Geo.Li)+Geo.Li; %width of panel in m

%-------------------------------------------------------------------------%

%% COMSOL PROBE INFORMATION
%-------------------------------------------------------------------------%
Probe.domain = 4; %radius of air domain in meters
Probe.radius = 3; %radius of arc in meters
Probe.theta_min = pi*(10/180); %arc starting angle
Probe.theta_max = pi*(170/180); %arc end angle
Probe.Resolution = 181;
Probe.theta_vector = linspace(Probe.theta_min,Probe.theta_max,Probe.Resolution);
Probe.Coordinates(1,:) = Probe.radius*cos(Probe.theta_vector); %Probe x coordinates
Probe.Coordinates(2,:) = Probe.radius*sin(Probe.theta_vector); %Probe y coordinates

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

calcDC = questdlg("calculate diffusion coefficient?");

switch calcDC
    case 'Yes'

    %% CALCULATE DIFFUSION COEFFICIENT
    %-------------------------------------------------------------------------   
    % QRD
    delta.COMSOL1 = calculateDiffusionCoefficient(Ps_1,Probe.theta_vector);
    %QRM
    delta.COMSOL2 = calculateDiffusionCoefficient(Ps_2,Probe.theta_vector);
    % Flat plane
    delta.flatnum = calculateDiffusionCoefficient(Psflatnum,Probe.theta_vector);
    
    run("QRD_TMM.m")
    figure()
    plot(Freq.Vector,delta.COMSOL1,"LineWidth",2); %plot diffusion coefficient
    hold on
    plot(Freq.Vector,delta.COMSOL2,"LineWidth",2)
    hold on
    plot(Freq.Vector,delta.QRD,"LineWidth",2)
    hold on
    plot(Freq.Vector,delta.flatnum,"LineWidth",2)
    hold on
    plot(Freq.Vector,delta.f,"LineWidth",2,"LineStyle","--")
    ylim([0, 1])
    legend(["N=5 QRD (numerical)","N=5 QRM (numerical)","N=5 QRD (TMM)","flat plane (numerical)","flat plane (TMM)"],...
        "Location","southeast")
    title(['diffusion coefficient - $r=$' num2str(Probe.radius) ' m'])
    xlabel("Hz")
    ylabel("$\delta$")
    %-------------------------------------------------------------------------%
    case 'No'
        return

end

%% Function for calculating the diffusion coefficient

function delta = calculateDiffusionCoefficient(Ps,theta)
    %% CALCULATE DIFFUSION COEFFICIENT
    %-------------------------------------------------------------------------%
    SI = abs(Ps).^2; %SPL is converted to sound intensity.
    n_d = length(theta);
    
    SIsum = sum(SI,2);
    SIsq = sum(SI.^2,2);
    
    delta = (SIsum.^2 - SIsq)./((n_d-1)*(SIsq)); %diffusion coefficient

end
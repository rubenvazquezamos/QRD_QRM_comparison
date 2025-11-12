%% FREQUENCY
%-------------------------------------------------------------------------%
Freq.f_min = 250;                                  % Minimum Freq of interest
Freq.f_max = 4000;                               % Maximum Freq of interest
Freq.df = 50;                                    % Freq discretization
Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);    % Freq vector
Freq.Nf = numel(Freq.Vector);                   % Number of frequencies

%% GEOMETRY
%-------------------------------------------------------------------------%
Geo.numberWells = 5; %prime number used to generate s_n
Geo.wellWidth = 7e-2;
Geo.betweenWells = 1e-2;
Geo.panelLength = Geo.numberWells*(Geo.wellWidth+Geo.betweenWells)...
    +Geo.betweenWells; %width of panel in m
Geo.stockThickness = 12e-3; %stock thickness
%-------------------------------------------------------------------------%

%% COMSOL FILE INFORMATION
%-------------------------------------------------------------------------%
File.Path = [pwd,filesep,'Models'];
File.Tag1 = 'Comsol_QRD5';
File.Tag2 = 'flat_plane';
File.Tag3 = 'optgeo.mphbin';
File.Extension = '.mph';

%% COMSOL PROBE INFORMATION
%-------------------------------------------------------------------------%
Probe.domain = 4; %radius of air domain in meters
Probe.radius = 3; %radius of arc in meters
Probe.theta_min =  0; %arc starting angle
Probe.theta_max = pi; %arc end angle
Probe.Resolution = 181;
Probe.theta_vector = linspace(Probe.theta_min,Probe.theta_max,Probe.Resolution);
Probe.Coordinates(1,:) = Probe.radius*cos(Probe.theta_vector); %Probe x coordinates
Probe.Coordinates(2,:) = Probe.radius*sin(Probe.theta_vector); %Probe y coordinates
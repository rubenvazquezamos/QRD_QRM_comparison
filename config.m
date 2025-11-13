%% FREQUENCY
%-------------------------------------------------------------------------%
Freq.f_min = 250;
Freq.f_max = 4000;
Freq.df = 50;                                    % Freq discretization
Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);
Freq.Nf = numel(Freq.Vector);                   % Number of frequencies
Freq.designFreq = 500; %QRD design frequency

%% GEOMETRY
%-------------------------------------------------------------------------%
Geo.numberWells = 5; %prime number used to generate s_n
Geo.wellWidth = 7e-2; %also unit cell width
Geo.betweenWells = 1e-2;
Geo.panelLength = Geo.numberWells*(Geo.wellWidth+Geo.betweenWells)...
    +Geo.betweenWells;
Geo.panelWidth = 3e-2;
Geo.stockThickness = 12e-3; %stock used for milling internal components
Geo.baseThickness = 8e-3; %stock used for base
Geo.lidThickness = 3e-3; %stock used for lid

%-------------------------------------------------------------------------%

%% ANGLES
% TMM: -pi/2 to pi/2
% COMSOL: 0 to pi

Angle.restriction = 0; %restricts extreme observer angles
Angle.theta_min = -pi/2;
Angle.theta_max = pi/2;
Angle.Resolution = 181;
Angle.theta_Vector = linspace(Angle.theta_min,Angle.theta_max,Angle.Resolution);
%% FREQUENCY
%-------------------------------------------------------------------------%
Freq.f_min = 250;
Freq.f_max = 4000;
Freq.df = 50;                                    % Freq discretization
Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);
Freq.Nf = numel(Freq.Vector);                   % Number of frequencies
Freq.designFreq = 500; %QRD design frequency
Freq.fangular = 2*pi*Freq.Vector; % angular frequency w
Freq.waveno = Freq.fangular./c_0; %wavenumber k

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
Geo.numxpts = lcm(Geo.numberWells,100); %number of x points
Geo.xpoints = linspace(-Geo.panelLength./2,Geo.panelLength./2,Geo.numxpts); %x points

%-------------------------------------------------------------------------%

%% ANGLES
% TMM: -pi/2 to pi/2
% COMSOL: 0 to pi

Angle.restriction = 0; %restricts extreme observer angles
Angle.min = -pi/2;
Angle.max = pi/2;
Angle.Resolution = 181;
Angle.Vector = linspace(Angle.min,Angle.max,Angle.Resolution);
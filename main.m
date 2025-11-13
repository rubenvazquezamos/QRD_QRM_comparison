% QRD QRM COMPARISON
% 
% This script optimises a geometry and later generates it in COMSOL
% to compare diffusion coefficients.
% % Scripts are Loosely Coupled, so variables are may be restated
% across scripts. config.m file ensures consistency of some parameters.
%
% HOW TO USE
% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.
% 

%% OPTIMISE GEOMETRY
LatexPreamble()
run("config.m")
run("optimiseGeometry.m")
optgeo.panel_thickness = Geo.panelDepth;
optgeo.stock_thickness = Geo.stockThickness;
save("optgeo.mat","optgeo") %save optimised geometry
clear all %information hiding

%% GENERATE COMSOL GEOMETRY
load("optgeo.mat")
run("config")
geometry = prepGeometryforCOMSOL(optgeo);
QRMgeometryCOMSOL = generateCOMSOLgeom(geometry,Geo.numberWells,Geo.stockThickness);
QRMgeometryCOMSOL.resetHist %compact model history
QRMgeometryCOMSOL.geom('geom1').export('optgeo.mphbin') %geometry file for import

%% OBTAIN SCATTERED PRESSURE
run("buildCOMSOLmodels.m")
run("QRD_TMM.m")
run("QRM_TMM.m")
run("flat_TMM.m")

%% DIFFUSION COEFFICIENT

model_names = {'QRD_TMM', 'QRD_COMSOL', 'QRM_TMM', 'QRM_COMSOL', 'flat_TMM', 'flat_COMSOL'};

for k = 1:numel(model_names)
    delta.(model_names{k}) = calculateDiffusionCoefficient(Ps.(model_names{k}));
    deltan.(model_names{k}) = normaliseDiffusionCoefficient(delta.(model_names{k}));
end

%% RESULTS
for k = 1:numel(model_names)
    plot(Freq.Vector,delta.(model_names{k}))
    hold on
end
hold off
%% ----               FUNCTIONS              ----

function geometry = prepGeometryforCOMSOL(optgeo)
    
    geometry = struct();
    geometry.neck.l_n = optgeo.neck.lengths;
    geometry.neck.w_n = optgeo.neck.widths;
    geometry.cavity.w_c = optgeo.cavity.widths;
    geometry.cavity.l_c = optgeo.cavity.lengths;
    geometry.slit.hh = optgeo.slit.widths;
    geometry.slit.a_y = optgeo.slit.lengths;
    geometry.L = optgeo.panel_thickness;
    geometry.e = optgeo.stock_thickness;

end

function [] = LatexPreamble(~) % Latex Preamble
    set(groot,'defaultAxesTickLabelInterpreter','latex');  
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultAxesFontSize',18)
    set(0,'DefaultFigureWindowStyle','docked');
    set(0,'defaultFigureColor',[1 1 1])
    path = convertCharsToStrings(fileparts(matlab.desktop.editor.getActiveFilename));
    cd(path)
    clear; clear global *; clc; warning off;
end


function Rbig = discretizeDiffuserReflectionCoefficient(N,x,R)
    % INPUTS
    %   N: number of wells and prime number
    %   x: vector of points along diffuser surface.
    %   R: matrix of reflection coefficient values.
    % 
    % OUTPUTS
    %   Rbig: matrix of spatially discretized reflection coefficient values.
            
        augmentation_index = length(x)./N;
        Rbig = repelem(R, 1, augmentation_index);    
end

function Ps = fftfraunhofer(theta,Rs,k,x)
    % ======================================
    %
    %   PS = FFTFRAUNHOFER(THETA,RS,K,X)
    %
    %   FFTFRAUNHOFER calculates the far field PS using the Fraunhofer integral
    %   of a surface with reflection RS(X) at wavenumber k and at angles
    %   theta
    %
    %   Noé Jiménez, Salford, October 2016
    %
    %=======================================
    na = length(x);
    dx = x(2)-x(1);
    Ps = zeros(size(theta));
    for ia=1:na
        Ps = Ps+Rs(ia).*exp(1i*k*x(ia)*sin(theta))*dx; 
    end

end

function delta = calculateDiffusionCoefficient(Ps,theta,ri)
    % INPUTS:
    %   Ps: [NxM] matrix of scattered pressure across N wavenumbers 
    %       and N observer angles
    %   theta: vector of observer angles
    %   ri: scalar restriction on extreme observer angles
    %      eg: ri = 2 gives angle range -88 to 88 deg.
    % OUTPUT:
    %   delta: vector of diffusion coefficient values
    %-------------------------------------------------------------------------%
   
    if ri ~= 0
      theta = Geo.Theta(ri:end-ri); %restrict extreme observer angles
      Ps = Ps(:,ri:end-ri);
    end

    SI = abs(Ps).^2; %SPL is converted to sound intensity.
    n_d = length(theta);
    
    SIsum = sum(SI,2);
    SIsq = sum(SI.^2,2);
    
    delta = (SIsum.^2 - SIsq)./((n_d-1)*(SIsq)); %diffusion coefficient

end

function delta_n = normaliseDiffusionCoefficient(delta, deltaf)
    %INPUTS
    % delta: diffusion coefficient of reflector
    % deltaf: diffusion coeffficient of flat plane
    %
    %OUTPUT
    % delta_n: normalised diffusion coefficient

    delta_n = (delta - deltaf)./(1-deltaf);
end
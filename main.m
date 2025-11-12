% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.

% This script optimises a geometry and later generates it in COMSOL
% to compare diffusion coefficients
% 
% Scripts are Loosely Coupled, so variables are may be restated
% across scripts.
% config.m file ensures consistency of some parameters across scripts

%% EXECUTE SCRIPTS
LatexPreamble()
run("config.m")
run("optimiseGeometry.m")

clear all %information hiding
load("optgeo.mat")

run("config")
geometry = structureGeometry(optgeo);
QRM_geometry = generateCOMSOLgeom(geometry,Geo.numberWells,Geo.stockThickness);
QRM_geometry.resetHist %compact model history
QRM_geometry.geom('geom1').export('optgeo.mphbin') %geometry file for import

run("buildCOMSOLmodels.m");

%% PLOTS



%% FUNCTIONS

function geometry = structureGeometry(optgeo)
    
    geometry = struct();
    geometry.neck.l_n = optgeo.neck.lengths;
    geometry.neck.w_n = optgeo.neck.widths;
    geometry.cavity.w_c = optgeo.cavity.widths;
    geometry.cavity.l_c = optgeo.cavity.lengths;
    geometry.slit.hh = optgeo.slit.widths;
    geometry.slit.a_y = optgeo.slit.lengths;
    geometry.L = optgeo.panel_thickness;

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

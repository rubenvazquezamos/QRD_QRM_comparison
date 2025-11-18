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
run("optimConfig.m")
rng("shuffle") %initialise random number generator
run("optimiseGeometry.m")
save("optgeo.mat","optgeo") %save optimised geometry
disp(initialguess)

return

clear all %information hiding

%% GENERATE COMSOL GEOMETRY
run("config")
load("optgeo.mat")
geoCOMSOL = prepCOMSOLgeom(optgeo,Geo);
QRMgeometryCOMSOL = generateCOMSOLgeom(geoCOMSOL);
QRMgeometryCOMSOL.resetHist %compact model history
QRMgeometryCOMSOL.geom('geom1').export('optgeo.mphbin') %geometry file for import

%% OBTAIN SCATTERED PRESSURE
run("buildCOMSOLmodels.m")
run("TMM_models.m")

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


%% ----               HELPER FUNCTIONS              ----

function geometry = prepCOMSOLgeom(optgeo,Geo)
    
    geometry = struct();

    geometry.neck.lengths = optgeo.neck.lengths;
    geometry.neck.widths = optgeo.neck.widths;
    geometry.cavity.widths = optgeo.cavity.widths;
    geometry.cavity.lengths = optgeo.cavity.lengths;
    geometry.slit.widths = optgeo.slit.widths;
    geometry.slit.lengths = optgeo.slit.lengths;

    % add extra required parameters
    geometry.stockThickness = Geo.stockThickness;
    geometry.panelLength = Geo.panelLength;
    geometry.panelWidth = Geo.panelWidth;
    geometry.numberWells = Geo.numberWells;
    geometry.baseThickness = Geo.baseThickness;
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
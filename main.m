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

clear all %information hiding

%% GENERATE COMSOL GEOMETRY
run("config")
load("optgeo.mat","optgeo")
geoCOMSOL = prepCOMSOLgeom(optgeo,Geo);
QRMgeometryCOMSOL = generateCOMSOLgeom(geoCOMSOL);
QRMgeometryCOMSOL.resetHist %compact model history
QRMgeometryCOMSOL.geom('geom1').export('optgeo.mphbin') %geometry file for import

%% SCATTERED PRESSURE
run("buildCOMSOLmodels.m")
run("TMM_models.m")
%% DIFFUSION COEFFICIENT
model_names = ["QRD TMM", "QRM TMM", "flat TMM","QRD numerical",...
    "QRM numerical","flat numerical"];

delta.QRD_TMM = calculateDiffusionCoefficient(Ps.QRD_TMM,Angle.Vector,5);
delta.QRM_TMM = calculateDiffusionCoefficient(Ps.QRM_TMM,Angle.Vector,5);
delta.flat_TMM = calculateDiffusionCoefficient(Ps.flat_TMM,Angle.Vector,5);
delta.QRD_COMSOL = calculateDiffusionCoefficient(Ps.QRD_COMSOL,Angle.Vector,5);
delta.QRM_COMSOL = calculateDiffusionCoefficient(Ps.QRM_COMSOL,Angle.Vector,5);
delta.flat_COMSOL = calculateDiffusionCoefficient(Ps.flat_COMSOL,Angle.Vector,5);

% 
% deltan.QRD_TMM = normaliseDiffusionCoefficient(delta.QRD_TMM,delta.flat_TMM);
% deltan.QRM_TMM = normaliseDiffusionCoefficient(delta.QRM_TMM,delta.flat_TMM);

%% RESULTS
figure()
set(groot,'DefaultLineLineWidth',2)
plot(Freq.Vector,delta.QRD_TMM)
hold on
plot(Freq.Vector,delta.QRM_TMM)
plot(Freq.Vector,delta.flat_TMM)
plot(Freq.Vector,delta.QRD_COMSOL)
plot(Freq.Vector,delta.QRM_COMSOL)
plot(Freq.Vector,delta.flat_COMSOL)
legend(model_names)

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
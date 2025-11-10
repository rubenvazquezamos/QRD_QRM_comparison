% 1 Run COMSOL to MATLAB LiveLink
% 2 Then run this script.

% This script optimises a geometry and later generates it in COMSOL
% to compare diffusion coefficients
% 
clear all
load("optgeo.mat")
% 
geometry = struct();
geometry.neck.l_n = optgeo.neck.lengths;
geometry.neck.w_n = optgeo.neck.widths;
geometry.cavity.w_c = optgeo.cavity.widths;
geometry.cavity.l_c = optgeo.cavity.lengths;
geometry.slit.hh = optgeo.slit.widths;
geometry.slit.a_y = optgeo.slit.lengths;
geometry.L = optgeo.panel_thickness;

N = 5;
e = 1e-3;

QRM_geometry = geo_generator(geometry,N,e); %generate geometry
QRM_geometry.resetHist %compact model history
QRM_geometry.geom('geom1').export('optgeo.mphbin') %save .mphbin geometry file

run("scatteredPressure_main");
% 
% %% ------ IMPORT METADIFFUSER GEOMETRY ------------------------------------
% model.component('comp1').geom('geom1').create('imp1', 'Import');
% model.component('comp1').geom('geom1').feature('imp1').set('type', 'native');
% model.component('comp1').geom('geom1').feature('imp1').set('filename', File.Tag3);
% model.component('comp1').geom('geom1').create('dif1', 'Difference');
% model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'c1'});
% model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'imp1'});
% %% ------------------------------------------------------------------------


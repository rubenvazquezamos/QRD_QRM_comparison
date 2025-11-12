clear all

%% COMSOL Acoustic Metadiffuser Geometry Generator
% This script calls the function unit_cell, which generates COMSOL geometry
% for an acoustic metadiffuser.
% Step 1: Run LiveLink
% Step 2: run this script

geometry = struct();
geometry.cavity.l_c = [10 20 30];
geometry.cavity.w_c = [20 20 20];
geometry.neck.l_n = [10 20 30];
geometry.neck.w_n = [5 10 15];
geometry.slit.hh = [10 20 30];
geometry.slit.a_y = [21 21 21];
geometry.L = 22;

geo_generator(geometry,3,1)
function [pressure] = QR_5(Freq,Geo,Probe,File,save_dlg)
%
% QR_5_test1.m
%
% Model exported on Oct 12 2025, 14:26 by COMSOL 6.0.0.318.

%-------------------------------------------------------------------------%
%% CALL COMSOL
disp(' -- Call COMSOL')
%-------------------------------------------------------------------------%

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');


model.component.create('comp1', true);


%-------------------------------------------------------------------------%
%% PARAMETERS
disp(' -- Sending parameters')
%-------------------------------------------------------------------------%

model.param.set('L', num2str(Geo.panelLength), 'Length of the diffuser');
model.param.set('Lw', num2str(Geo.wellWidth), 'Width of one well');
model.param.set('Li', num2str(Geo.betweenWells), 'Width in between wells');
model.param.set('H', '0.3[m]', 'Height of the diffuser');
model.param.set('d1', '0[m]', 'Depth of well 1');
model.param.set('d2', '0.068[m]', 'Depth of well 2');
model.param.set('d3', '0.272[m]', 'Depth of well 3');
model.param.set('d4', '0.272[m]', 'Depth of well 4');
model.param.set('d5', '0.068[m]', 'Depth of well 5');
model.param.set('r_air', num2str(Probe.domain), 'Radius of the air domain (for single diffuser)');
model.param.set('Hpml', '0.2[m]', 'Thickness of the PML');
model.param.group.create('par2');
model.param('par2').set('c0', '343[m/s]', 'Speed of sound');
model.param('par2').set('rho0', '1.225[kg/m^3]', 'Density');
model.param('par2').set('f0', '125[Hz]', 'Frequency parameter');
model.param('par2').set('f_min', num2str(Freq.f_min));
model.param('par2').set('f_max', num2str(Freq.f_max));
model.param('par2').set('df', num2str(Freq.df));
model.param('par2').set('phi', '0 [rad]', 'phase angle');
model.param.label('Geometry');
model.param('par2').label('Physics');

%-------------------------------------------------------------------------%
%% COMSOL PHYSICAL CONSTANTS
disp(' -- Set global constants')
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%% VARIABLES
disp(' -- Sending variables')
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%% GEOMETRY
disp(' -- Build geometry')
%-------------------------------------------------------------------------%

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('layername', {'PML'});
model.component('comp1').geom('geom1').feature('c1').setIndex('layer', '0.75', 0);
model.component('comp1').geom('geom1').feature('c1').set('r', 'r_air');
model.component('comp1').geom('geom1').feature('c1').set('angle', 180);
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').active(false);
model.component('comp1').geom('geom1').feature('r1').set('pos', {'-L/2+Li' '-d1'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'Lw' 'd1'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'-L/2+2*Li+Lw' '-d2'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'Lw' 'd2'});
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'-L/2+3*Li+2*Lw' '-d3'});
model.component('comp1').geom('geom1').feature('r3').set('size', {'Lw' 'd3'});
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'-L/2+4*Li+3*Lw' '-d4'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'Lw' 'd4'});
model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').set('pos', {'-L/2+5*Li+4*Lw' '-d5'});
model.component('comp1').geom('geom1').feature('r5').set('size', {'Lw' 'd5'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').set('p', {'-L/2+Li' '0'});
model.component('comp1').geom('geom1').nodeGroup.create('grp1');
model.component('comp1').geom('geom1').nodeGroup('grp1').label('Wells');
model.component('comp1').geom('geom1').nodeGroup('grp1').placeAfter('c1');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r1');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r2');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r3');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r4');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r5');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');


%%% EXTRA CODE: visualise the geometry before going further with
% approval input with the [Enter] key.
fig = figure();
clf(fig)
set(fig,'Renderer','opengl');
mphgeom(model,'geom1','facealpha',0.5);
box on;
% kk = input('Is the geometry valid?');clear kk;


%-------------------------------------------------------------------------%
%% PHYSICS
disp(' -- Implementing physics')
%-------------------------------------------------------------------------%
model.component('comp1').physics.create('acpr', 'PressureAcoustics', 'geom1');
model.component('comp1').physics('acpr').create('bpf1', 'BackgroundPressureField', 2);
model.component('comp1').physics('acpr').feature('bpf1').selection.set([1 2 3 4 5 6 7]);
model.component('comp1').physics('acpr').create('imp1', 'Impedance', 1);
model.component('comp1').physics('acpr').feature('imp1').selection.set([2 25]);


model.component('comp1').physics('acpr').feature('bpf1').set('p', 'p_inc');
model.component('comp1').physics('acpr').feature('bpf1').set('dir', [0; -1; 0]);
model.component('comp1').physics('acpr').feature('bpf1').set('phi', 'phi');
model.component('comp1').physics('acpr').feature('bpf1').set('pamp', 1);
model.component('comp1').physics('acpr').feature('bpf1').set('c_mat', 'from_mat');
%-------------------------------------------------------------------------%
%% PML
disp(' -- Set PML')
%-------------------------------------------------------------------------%

model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.set([1 5]);
model.component('comp1').coordSystem('pml1').set('ScalingType', 'Cylindrical');


%-------------------------------------------------------------------------%
%% MATERIAL
disp(' -- Define materials')
%-------------------------------------------------------------------------%

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup.create('NonlinearModel', 'Nonlinear model');
model.component('comp1').material('mat1').propertyGroup.create('idealGas', 'Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func.create('Cp', 'Piecewise');


model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '293.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
model.component('comp1').material('mat1').materialType('nonSolid');

%-------------------------------------------------------------------------%
%% STUDY
disp(' -- Creating study')
%-------------------------------------------------------------------------%

model.study.create('std1');
model.study('std1').create('freq', 'Frequency');
model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.study('std1').feature('freq').set('plist', 'range(f_min,df,f_max)');

%-------------------------------------------------------------------------%
%% MESH
disp(' -- Build mesh')
%-------------------------------------------------------------------------%

model.component('comp1').mesh.create('mesh1');

%-------------------------------------------------------------------------%
%% SOLVE
disp(' -- Solving')
%-------------------------------------------------------------------------%
%%% EXTRA CODE: opens a window to show the solver's progress
ModelUtil.showProgress(true);

model.sol('sol1').runAll;
%-------------------------------------------------------------------------%
%% PLOT
disp(' -- Generate plots')
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%% GET DATA/PROBING
disp(' -- Probing')
%-------------------------------------------------------------------------%
pressure = mphinterp(model,'acpr.p_s','coord',Probe.Coordinates);

%-------------------------------------------------------------------------%
%% SAVE
disp(' -- Save model')
%-------------------------------------------------------------------------%
% 

if save_dlg == true 
    save_dlg = questdlg("clear solution and mesh data?");
    switch save_dlg
        case 'Yes'
        model.sol('sol1').clearSolutionData; %Clear solution data
        model.mesh.clearMeshes; %Clear meshes
        case 'No' 
    end
end

mphsave(model,[File.Path,filesep,File.Tag1,File.Extension]);

disp(' -- DONE')

end


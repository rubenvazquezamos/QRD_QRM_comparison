function out = model
%
% QRM_test.m
%
% Model exported on Nov 10 2025, 12:04 by COMSOL 6.0.0.318.

%-------------------------------------------------------------------------%
%% CALL COMSOL
disp(' -- Call COMSOL')
%-------------------------------------------------------------------------%

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\User\Desktop\comparison 2');

model.label('QRM.mph');

%-------------------------------------------------------------------------%
%% PARAMETERS
disp(' -- Sending parameters')
%-------------------------------------------------------------------------%

model.param.set('L', '0.5[m]', 'Width of the diffuser');
model.param.set('r_air', '3', 'Radius of the air domain (for single diffuser)');
model.param.set('r0', '100[m]', 'Evaluation distance');
model.param.set('Hair', '1[m]', 'Height of the air domain');
model.param.set('Hpml', '0.2[m]', 'Thickness of the PML');
model.param.group.create('par2');
model.param('par2').set('c0', '343[m/s]', 'Speed of sound');
model.param('par2').set('rho0', '1.225[kg/m^3]', 'Density');
model.param('par2').set('f_min', '125');
model.param('par2').set('f_max', '250');
model.param('par2').set('df', '5');
model.param('par2').set('phi', '0 [rad]', 'phase angle');
model.param.label('Geometry');
model.param('par2').label('Physics');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('layername', {'PML'});
model.component('comp1').geom('geom1').feature('c1').setIndex('layer', '0.75', 0);
model.component('comp1').geom('geom1').feature('c1').set('r', 'r_air');
model.component('comp1').geom('geom1').feature('c1').set('angle', 180);
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').active(false);
model.component('comp1').geom('geom1').feature('uni1').set('intbnd', false);
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'c1'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').set('p', {'-L/2+Li' '0'});
%% ------ IMPORT METADIFFUSER GEOMETRY ------------------------------------
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('type', 'native');
model.component('comp1').geom('geom1').feature('imp1').set('filename', 'C:\Users\User\Desktop\comparison 2\optgeo.mphbin');
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'c1'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'imp1'});
%% ------------------------------------------------------------------------
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('p_inc', '1[Pa]*exp(i*2*pi*freq/c0*(sin(theta0)*x-cos(theta0)*y))', 'Incident plane wave');
model.component('comp1').variable('var1').set('p_inf', '1[Pa]*exp(i*2*pi*freq/c0*(sin(theta0)*x+cos(theta0)*y))', 'Reflected plane wave from infinite baffle');
model.component('comp1').variable('var1').set('p_scat', 'p_inf+acpr.p_s', 'Scattered pressure');
model.component('comp1').variable('var1').set('p_scat_ext', 'pext(x,y)', 'Scattered pressure in the exterior field');

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

model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.set([1 3]);

model.component('comp1').physics.create('acpr', 'PressureAcoustics', 'geom1');
model.component('comp1').physics('acpr').create('bpf1', 'BackgroundPressureField', 2);
model.component('comp1').physics('acpr').feature('bpf1').selection.set([1 2 3]);
model.component('comp1').physics('acpr').create('pmb1', 'PerfectlyMatchedBoundary', 1);
model.component('comp1').physics('acpr').feature('pmb1').selection.set([70 73]);
model.component('comp1').physics('acpr').create('pmb2', 'PerfectlyMatchedBoundary', 1);
model.component('comp1').physics('acpr').feature('pmb2').selection.set([1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69]);
model.component('comp1').physics('acpr').create('imp1', 'Impedance', 1);
model.component('comp1').physics('acpr').feature('imp1').selection.set([2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68]);

model.component('comp1').mesh('mesh1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').create('se1', 'SizeExpression');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('bl1', 'BndLayer');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size1').selection.set([1 2 3]);
model.component('comp1').mesh('mesh1').feature('se1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('se1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([1 2 3]);
model.component('comp1').mesh('mesh1').feature('bl1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('bl1').selection.set([1 2 3]);
model.component('comp1').mesh('mesh1').feature('bl1').create('blp1', 'BndLayerProp');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp1').selection.set([70 73]);
model.component('comp1').mesh('mesh1').feature('map1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('map1').selection.set([1 3]);
model.component('comp1').mesh('mesh1').feature('map1').create('dis1', 'Distribution');

model.component('comp1').view('view1').axis.set('xmin', -3.354295492172241);
model.component('comp1').view('view1').axis.set('xmax', 3.1817524433135986);
model.component('comp1').view('view1').axis.set('ymin', -1.1280899047851562);
model.component('comp1').view('view1').axis.set('ymax', 3.262068271636963);

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

model.component('comp1').coordSystem('pml1').set('ScalingType', 'Cylindrical');

model.component('comp1').physics('acpr').feature('bpf1').set('p', 'p_inc');
model.component('comp1').physics('acpr').feature('bpf1').set('dir', [0; -1; 0]);
model.component('comp1').physics('acpr').feature('bpf1').set('phi', 'phi');
model.component('comp1').physics('acpr').feature('bpf1').set('pamp', 1);
model.component('comp1').physics('acpr').feature('bpf1').set('c_mat', 'from_mat');
model.component('comp1').physics('acpr').feature('pmb1').active(false);
model.component('comp1').physics('acpr').feature('pmb2').set('directionType', 'normal');
model.component('comp1').physics('acpr').feature('pmb2').active(false);

model.component('comp1').mesh('mesh1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size1').set('hmax', 'c0/f_max/5');
model.component('comp1').mesh('mesh1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size1').set('hmin', 'c0/f_max/10');
model.component('comp1').mesh('mesh1').feature('size1').set('hminactive', false);

model.sol.create('sol1');

model.study.create('std1');
model.study('std1').create('freq', 'Frequency');
model.study('std1').feature('freq').set('useadvanceddisable', true);

model.component('comp1').mesh('mesh1').feature('se1').set('evaltype', 'initialexpression');
model.component('comp1').mesh('mesh1').feature('se1').set('sizeexpr', 'subst(real(acpr.c_c),acpr.freq,freqmax)/freqmax/5');
model.component('comp1').mesh('mesh1').feature('se1').set('studystep', 'std1/freq');
model.component('comp1').mesh('mesh1').feature('bl1').active(false);
model.component('comp1').mesh('mesh1').feature('bl1').set('sharpcorners', 'trim');
model.component('comp1').mesh('mesh1').feature('bl1').set('smoothtransition', false);
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp1').set('blnlayers', 3);
model.component('comp1').mesh('mesh1').feature('map1').active(false);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').set('numelem', 3);

model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('pDef', 'Parametric');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', 'abs(acpr.p_s)');
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').feature('surf1').set('expr', 'acpr.Lp_t');

model.study('std1').feature('freq').set('plist', 250);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Frequency Domain');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'250[Hz]'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('s1').feature('pDef').label('Parametric 2');
model.sol('sol1').feature('s1').feature('pDef').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('pDef').set('plistarr', [250]);
model.sol('sol1').feature('s1').feature('pDef').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('pDef').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('pDef').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('pDef').set('uselsqdata', false);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', [250]);
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').runAll;

model.result('pg1').label('Acoustic Pressure (acpr)');
model.result('pg1').set('showlegendsunit', true);
model.result('pg1').feature('surf1').set('colortable', 'ThermalClassic');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label('Sound Pressure Level (acpr)');
model.result('pg2').set('showlegendsunit', true);
model.result('pg2').feature('surf1').set('resolution', 'normal');

out = model;

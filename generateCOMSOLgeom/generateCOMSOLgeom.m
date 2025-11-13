function out = generateCOMSOLgeom(geometry)
%
% This function generates an Acoustic Metadiffuser geometry with
% a specified number of unit cells.
% cells.
%Inputs:
% geometry - strut() with fields
%   slit, cavity, neck, (with subfields 'widths' and 'lengths'.) 
%   stockThickness
%   lidThickness
%   baseThickness
%   panelLength
%   panelWidth
%   numberWells  
%
% Output
%
% COMSOL model of Metadiffuser geometry
%

%% MATLAB variables
e_stock = geometry.stockThickness;
e_base = geometry.baseThickness;

N = geometry.numberWells;

%rectangle labels

for r_index = 1:3*N
    r_label(r_index) = "r"+ num2str(r_index);
end

offset = geometry.stockThickness; %initialise cavity x position
D = 0; %initialise panel length
L = geometry.panelWidth + geometry.baseThickness; %recall w_c = L - lidThickness
%% COMSOL MODEL

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

% model.modelPath('C:\Users\User\Dropbox\PhD\Projects\COMSOL_geo_parametric');

model.label('metadiffuser.mph');

model.param.group.create('par2');
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 2);


for ii = 1:N %rectangle loop
    
    l_c = geometry.cavity.lengths(ii);
    w_c = geometry.cavity.widths(ii);
    l_n = geometry.neck.lengths(ii);
    w_n = geometry.neck.widths(ii);
    hh = geometry.slit.widths(ii);
    a_y = geometry.slit.lengths(ii);
    

    index = 3.*(ii-1) + 1;

    %dimensions
    model.param.set(['w_c', num2str(ii)], num2str(w_c,'%.10f') , 'width cavity');
    model.param.set(['l_c', num2str(ii)], num2str(l_c,'%.10f'), 'length cavity');
    model.param.set(['w_n', num2str(ii)], num2str(w_n,'%.10f'), 'width neck');
    model.param.set(['l_n', num2str(ii)], num2str(l_n,'%.10f'), 'length neck');
    model.param.set(['a_y', num2str(ii)], num2str(a_y,'%.10f'), 'slit width');
    model.param.set(['hh', num2str(ii)], num2str(hh,'%.10f'), 'slit length');
    model.param.set('e_stock', num2str(e_stock,'%.10f'), 'stock thickness');
    model.param.set('e_base', num2str(e_base,'%.10f'), 'base thickness');
    model.param.set('L', [num2str(L)], 'panel width');
    model.param.label('dimensions');
    
    %% Definitions
    %Cavity
    model.component('comp1').geom('geom1').create( num2str(r_label(index)) , 'Rectangle');
    model.component('comp1').geom('geom1').feature( num2str(r_label(index)) ).label(['cavity', num2str(ii)]);
    model.component('comp1').geom('geom1').feature( num2str(r_label(index)) ).set('pos', {num2str(offset,'%.10f') 'e_base'});
    model.component('comp1').geom('geom1').feature( num2str(r_label(index)) ).set('size', {['l_c',num2str(ii)], ['w_c',num2str(ii)]});
    
    offset = ['l_c',num2str(ii), ' + ', num2str(offset,'%.10f')]; %augment offset
    
    mphgeom(model,'geom1','facealpha',0.5);
    %Neck
    model.component('comp1').geom('geom1').create( num2str(r_label(index+1)) , 'Rectangle');
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+1)) ).label(['neck', num2str(ii)]);
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+1)) ).set( 'pos', {...
            [num2str(offset,'%.10f')] ['(w_c', num2str(ii),'-w_n',num2str(ii),')/2+e_base']...
            } );
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+1)) ).set('size', {['l_n',num2str(ii)] ['w_n',num2str(ii)]});
    
    offset = ['l_n',num2str(ii),' + ', num2str(offset,'%.10f')]; %augment offset
    
    mphgeom(model,'geom1','facealpha',0.5);
    %Slit
    model.component('comp1').geom('geom1').create( num2str(r_label(index+2)) , 'Rectangle');
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+2)) ).label(['slit', num2str(ii)]);
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+2)) ).set('pos', {num2str(offset,'%.10f') 'e_base'});
    model.component('comp1').geom('geom1').feature( num2str(r_label(index+2)) ).set('size', {['hh',num2str(ii)] ['a_y',num2str(ii)]});
    
    offset = ['e_stock + hh',num2str(ii), '+', num2str(offset,'%.10f')]; %augment offset
    
    %enclosing box
    D = ['l_c',num2str(ii),'+ l_n',num2str(ii),'+ hh',num2str(ii), '+ e_stock +', num2str(D)];  %augment panel length
    model.param.set('D', num2str(D,'%.10f'), 'panel length'); %set panel size
    mphgeom(model,'geom1','facealpha',0.5);
end
D = ['e_stock +', num2str(D)];  %final agumentation
model.param.set('D', num2str(D,'%.10f'), 'panel length'); %set panel size

mphgeom(model,'geom1','facealpha',0.5);
%% Build enclosing rectangle
r_label(index+3) = ['r', num2str(r_index+1)];
%Panel
model.component('comp1').geom('geom1').create(num2str(r_label(index+3)), 'Rectangle');
model.component('comp1').geom('geom1').feature( num2str(r_label(index+3)) ).set('size', {'D' 'L'});
model.component('comp1').geom('geom1').feature( num2str(r_label(index+3)) ).set('pos', {'0' '0'});

inner_geometry = convertStringsToChars( r_label(1:(end-1)) );
outer_geometry = convertStringsToChars( r_label(end) );

model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').set('intbnd', false);
model.component('comp1').geom('geom1').feature('uni1').selection('input').set(inner_geometry);

mphgeom(model,'geom1','facealpha',0.5);
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').selection('input').set(outer_geometry);

model.component('comp1').geom('geom1').feature('dif1').selection('input2').set('uni1');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

% center the panel at x = 0
model.component('comp1').geom('geom1').create('mov1', 'Move');
model.component('comp1').geom('geom1').feature('mov1').set('displx', '-D/2');
model.component('comp1').geom('geom1').feature('mov1').selection('input').set({'dif1'});

mphgeom(model,'geom1','facealpha',0.5);

out = model;

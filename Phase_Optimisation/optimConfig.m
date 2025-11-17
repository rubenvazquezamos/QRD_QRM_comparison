%% OPTIMISATION CONFIGURATION
% min values should be entered as -min into matrix
% max values are positive

options = optimoptions('fmincon','Display','iter','PlotFcn',{@optimplotstepsize...
,@optimplotfval,@optimplotx},'MaxFunctionEvaluations',...
500*length(initialguess(:)),'MaxIterations',500*length(initialguess(:)),...
'FunctionTolerance',1e-7,'StepTolerance',1e-14,'Algorithm','interior-point');



a_x = Geo.panelLength./Geo.numberWells; %unit cell length
L = Geo.panelWidth;
e_l = Geo.lidThickness;
e_s = Geo.stockThickness;

minw_n = L - 2*e_s;
maxw_n = L - (2/3)*e_s;
minl_n = 15e-3;
minl_c = 1e-3;
minh = 1e-3;

% Inequality constraints
% [w_n,l_n,w_c,l_c,h,a_y]
A = [
      1  0 -1  0  0  0;
     -1  0  0  0  0  0;
      1  0  0  0  0  0;
      0 -1  0  0  0  0;
      0  0  0 -1  0  0;
      0  0  0  0 -1  0;
     ];

b = [
    0;
    -minw_n;
    maxw_n;
    -minl_n;
    -minl_c;
    -minh;
    ];

% Equality constraints
% [w_n,l_n,w_c,l_c,h,a_y]
Aeq = [
        0 0 0 0 0 1;
        0 0 1 0 0 0
        0 1  0 1 1  0
                    ];

beq = [L; L-e_l, ; a_x];
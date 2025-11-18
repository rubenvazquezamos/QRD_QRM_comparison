%This code matches the phase of the reflection coefficient of a single AMD
%well at a single frequency to that of a QRD.
% Based on Phase_Optimisation v 2.1.2
% needs config.m to run
% 
% functions defined locally due to need for global variables

global targetphase n optim_f df stinson_constants

stinson_constants = generateStinsonConstants();

% =========================================================================

%% Obtain target phase
rho_0 = stinson_constants.density;
c_0 = stinson_constants.sound_speed;
n = Geo.numberWells;
df = 1; %frequency step
optim_f = 2000; %optimisation frequency
[QRD_depths, s_n, Rgoal] = QRDReflectionCoefficient(Geo.numberWells,...
    Geo.wellWidth,Freq.designFreq,optim_f,optim_f,df,rho_0,c_0);
targetphase = angle(Rgoal);

%---- Check phase profile ----
figure(1)
plotPhaseStep(n,angle(Rgoal))
title("$\arg(R_{target})$")

answer = questdlg('Desired target phase profile?');
switch answer
    case 'Yes' %script continues
    case 'No'
     error('Desired phase profile was not generated');
 end

% ========================================================================
%% Optimisation setup
run('optimConfig.m')

% ---- Build full constraint matrices ----
inequality_coefficients = [];
equality_coefficients = [];

for i = 1:n
    % Inequality constraints for well i
    row_ineq = zeros(size(A,1), 6*n); %define well matrix
    row_ineq(:, (i-1)*6 + 1 : i*6) = A; %populate with constaints
    inequality_coefficients = [inequality_coefficients; row_ineq]; %concatenate with other wells
    
    % Equality constraints for well i
    row_eq = zeros(size(Aeq,1), 6*n);
    row_eq(:, (i-1)*6 + 1 : i*6) = Aeq;
    equality_coefficients = [equality_coefficients; row_eq];
end

inequality_constants = repmat(b, n, 1);
equality_constants = repmat(beq, n, 1);

% Define bounds
upper_bounds = 50*ones(6*n,1)*1e-3;
lower_bounds = upper_bounds./1000;

% Initial guess that satisfies constraints
initialguess = (lower_bounds + upper_bounds)./2;

% ========================================================================
%% Get user input

eqconstraints = questdlg("enable equality constraints?");
switch eqconstraints
    case 'Yes'
    case 'No'
    equality_coefficients = [];
    equality_constants = [];
end

ineqconstraints = questdlg("enable inequality constraints?");
switch ineqconstraints
    case 'Yes'
    case 'No'
    inequality_coefficients = [];
    inequality_constants = [];
end

genetic = questdlg("run genetic algorithm to find initial guess?");
switch genetic
    case 'Yes'
        % positions: [w_n,l_n,w_c,l_c,h,a_y]
        initialguess = ga(@objectiveFunction,30,[],...
            [],equality_coefficients,...
            equality_constants,lower_bounds,upper_bounds);
        disp(initialguess);
    case 'No'

        randomise = questdlg("randomise initial guess?");
        switch randomise
        case 'Yes'
            initialguess = abs(rand(size(initialguess),"double"));
            initialguess = initialguess.*1e-3; %scale initial guess       
        case 'No'

            force = questdlg("force ballpark initial guess?");
            switch force
            case 'Yes'
                initialguess = cell2mat(flattenStruct2Cell(load('optgeo1.mat')));
                initialguess = initialguess';
                initialguess = initialguess(:); %ballpark initial guess
             case 'No'
            
            end

        end
end




% =======================================================================
%% Call fmincon to perform Optimisation

% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

[G_opt, fval] =fmincon(@objectiveFunction,initialguess,inequality_coefficients...
    ,inequality_constants,equality_coefficients,equality_constants,...
    lower_bounds,upper_bounds,[],options);

optgeo = unpackgeometry(G_opt,n);

% =======================================================================
%% Results
%---- Check optimised phase ----
figure(1)
hold on
phasecheck = angle(MetadiffuserReflectionCoefficient(n,...
    optgeo,optim_f,optim_f,df));
plotPhaseStep(n,phasecheck);
legend('goal','optimised')
title("$\arg(R_{optimised})$");
hold off

answer = questdlg('Desired optimised phase profile?');
switch answer
    case 'Yes' %script continues
    case 'No'
     return %script terminates if the desired phase profile is not generated.
end

%---- Check Geometry ----
table = geotable(optgeo);
disp(table)

figure()
spaceBetweenCellsinGraphic = 10e-3;
geomchecker(optgeo,n,spaceBetweenCellsinGraphic);
answer = questdlg('Feasible geometry?');

switch answer
    case 'Yes'
    case 'No'
        disp('geometry deemed infeasible by user')
        return
end

% ========================================================================
%% ----          OBJECTIVE FUNCTION           ----
function error = objectiveFunction(X)  %(X is geometry optimisation variable)
    global targetphase n optim_f df
    % x contains geometry parameters
    % Calculate QRD and metadiffuser phases
    %Unpack geometry structure
    geometry = unpackgeometry(X,n);
    
    metaR = MetadiffuserReflectionCoefficient(n,...
        geometry,optim_f,optim_f,df);
    metaphase = angle(metaR);
    % Objective is the absolute error between phases
    error = sum((metaphase - targetphase).^2);
end

% =========================================================================
%%               HELPER FUNCTIONS

function plotPhaseStep(n, phi)
    % inputs
    % phi = angle(R) is the phase in radians
    % n is the number of wells
    
    x = 0:n; %space vector 
    phi_ext = [phi, phi(end)];
    % Plot
    ticks = x + 0.5;
    stairs(x, phi_ext, 'LineWidth', 1.5);
    xticks(ticks)
    xticklabels(x+1)
    xlim([0,x(end)])

    xlabel('well number');
    ylabel('Phase [rad]');

    grid on;
end


function geometry = unpackgeometry(X,n)
%the function unpacks a geometry vector X into n wells.
% positions: [w_n,l_n,w_c,l_c,h,a_y]
    geometry = struct();
    ind = (2:6:(6*n-1)); %index for vectorised loop
    geometry.neck.widths = X(ind-1); %w_n
    geometry.neck.lengths = X(ind); % l_n
    geometry.cavity.widths = X(1+ind); % w_c
    geometry.cavity.lengths = X(2+ind); % l_c
    geometry.slit.widths = X(3+ind); % h
    geometry.slit.lengths = X(4+ind); % a_y 
end


function geomchecker(geometry,n, e)
    %creates graph of geometry for user validation
    % positions: [w_n,l_n,w_c,l_c,h,a_y] 
    
    figure()
    hold on
    % Variable initialisation
    offset = 0; %cumulative offset
    cavx = 0; %starting position
    
    for ii = 1:n
    
        w_n = geometry.neck.widths(ii); %y_dim
        l_n = geometry.neck.lengths(ii); %x_dim
        w_c = geometry.cavity.widths(ii); %y_dim
        l_c = geometry.cavity.lengths(ii); %x_dim
        h = geometry.slit.widths(ii); %x_dim
        a_y = geometry.slit.lengths(ii); %y_dim
        
        cavx = offset;
        cavy = 0;
        offset = offset+l_c;
        neckx = offset;
        necky = (w_c-w_n)./2;
        offset = offset+l_n;
        slitx = offset;
        slity = 0;
        offset = offset+h+e;
       
        rectangle('Position',[cavx cavy l_c w_c],'FaceColor',"#d4d4d0",'EdgeColor',[0 0 0]); %cavity
        rectangle('Position',[neckx necky l_n w_n],'FaceColor',"#d4d4d0",'EdgeColor',[0 0 0]); %neck
        rectangle('Position',[slitx slity h a_y],'FaceColor',"#d4d4d0",'EdgeColor',[0 0 0]); %slit
    end
        axis equal
        drawnow; % Force immediate plot update
        hold off
end

function table = geotable(optimisedgeometry)
    %outputs table of geometry values
    %geometry is a structure
    table = struct();

    table.w_n = optimisedgeometry.neck.widths';
    table.l_n = optimisedgeometry.neck.lengths';
    table.w_c = optimisedgeometry.cavity.widths';
    table.l_c = optimisedgeometry.cavity.lengths';
    table.h = optimisedgeometry.slit.widths';
    table.a_y = optimisedgeometry.slit.lengths';

    table = struct2table(table);

end

%==========================================================================
%%                      CALCULATION FUNCTIONS

function R = MetadiffuserReflectionCoefficient(n, geometry, f_min,f_max,df)
    % Function to calculate reflection coefficient
   
    global stinson_constants
    M = 1; %number of resonators per slit 
    f_v = f_min:df:f_max;
    w = 2*pi*f_v;
       
    %% Constants
    rho_0 = stinson_constants.density; %kg/m^3 air density
    B_0 = stinson_constants.bulk_modulus; %Pa adiabatic bulk modulus
    % Derived constantss
    Pr = stinson_constants.prandtl_no;
    eta = stinson_constants.kinematic_viscosity; %kinematic viscosity coeff.
    gam = stinson_constants.specific_heat_ratio; %ratio of specific heats

    %% Complex sound speed parameters
    G_rho = sqrt(1i*w*rho_0 / eta);
    G_k = sqrt(1i*w*rho_0*Pr / eta);
    
    for indn = 1:n % Loop over wells
        % Extract geometry for current well
        w_n = geometry.neck.widths(indn); %y_dim
        l_n = geometry.neck.lengths(indn); %x_dim
        w_c = geometry.cavity.widths(indn); %y_dim
        l_c = geometry.cavity.lengths(indn); %x_dim
        h = geometry.slit.widths(indn); %x_dim
        a_y = geometry.slit.lengths(indn); %y_dim
                
        S_0 = l_c + l_n + h;
        S_n = w_n;
        S_c = w_c;
        S_s = h;
        
        % Half-widths for calculations
        h_n = w_n / 2;
        h_c = w_c / 2;
        h_s = h / 2;
        
        % Calculate complex sound properties for each frequency
        %Bulk modulus, density and sound speed (vectorised)
        B_n = B_0 * (1 + (gam-1) * (tanh(h_n * G_k)) ./ (h_n * G_k)).^(-1);
        rho_n = rho_0 * (1 - (tanh(h_n * G_rho)) ./ (h_n * G_rho)).^(-1);
        c_n = sqrt(B_n./rho_n);
        
        B_c = B_0 * (1 + (gam-1) * (tanh(h_c * G_k)) ./ (h_c * G_k)).^(-1);
        rho_c = rho_0 * (1 - (tanh(h_c * G_rho)) ./ (h_c * G_rho)).^(-1);
        c_c = sqrt(B_c./rho_c);
        
        B_s = B_0 * (1 + (gam-1) * (tanh(h_s * G_k)) ./ (h_s * G_k)).^(-1);
        rho_s = rho_0 * (1 - (tanh(h_s * G_rho)) ./ (h_s * G_rho)).^(-1);
        c_s = sqrt(B_s./rho_s);
        
        % Calculate sound speeds
        c_n = sqrt(B_n./rho_n);
        c_c = sqrt(B_c./rho_c);
        c_s = sqrt(B_s./rho_s);
        c_0 = sqrt(B_0/rho_0);
        
        % Calculate impedances
        Z0 = (rho_0 * c_0) / S_0; % Impedance of free air
        Z_n = (rho_n .* c_n) ./ S_n; % Impedance of neck
        Z_c = (rho_c .* c_c) ./ S_c; % Impedance of cavity
        Z_s = (rho_s .* c_s) ./ S_s; % Impedance of waveguide
        
        % Calculate wavenumbers
        k_n = w ./ c_n;
        k_c = w ./ c_c;
        k_s = w ./ c_s;
        
        % Length corrections
        delta_l1 = 0.82 * (1 - 1.35 * (h_n / h_c) + 0.31 * (h_n / h_c)^3) * h_n;
        delta_l2 = 0.82 * (1 - 0.235 * (h_n / h_s) - 1.32 * (h_n / h_s)^2 + 1.54 * (h_n / h_s)^3 - 0.86 * (h_n / h_s)^4) * h_n;
        deltal = delta_l1 + delta_l2;
        
        % Radiation correction
        d = l_c + l_n + h;
        phi_t = h / d;
        m = 1:1:20; % Summation index
        delta_slit = h * phi_t * sum((sin(m*pi*phi_t).^2) ./ ((m*pi*phi_t).^3));
        
        % Slit impedance
        Z_slit = (1i * w * delta_slit * rho_0) / (phi_t * S_0);
        
        % Calculate Helmholtz resonator impedance
        Z_HR = 1i * Z_n .* (cos(k_n * l_n) .* cos(k_c * l_c) - (k_n * deltal .* Z_n ./ Z_c) .* cos(k_n * l_n) .* sin(k_c * l_c) - (Z_n ./ Z_c) .* sin(k_n * l_n) .* sin(k_c * l_c)) ...
             ./ (sin(k_n * l_n) .* cos(k_c * l_c) - (k_n * deltal .* Z_n ./ Z_c) .* sin(k_n * l_n) .* sin(k_c * l_c) + (Z_n ./ Z_c) .* cos(k_n * l_n) .* sin(k_c * l_c));
        
        % Waveguide transfer matrix elements
        T11 = cos(k_s * a_y/2);
        T12 = 1i * Z_s .* sin(k_s * a_y/2);
        T21 = 1i * sin(k_s * a_y/2) ./ Z_s;
        T22 = cos(k_s * a_y/2);
        
        % Calculate transfer matrices for each frequency
        for ifreq = 1:length(f_v)
            % Waveguide transfer matrix
            M_s = [T11(ifreq), T12(ifreq); T21(ifreq), T22(ifreq)];
            
            % HR transfer matrix
            M_HR = [1, 0; 1/Z_HR(ifreq), 1];
            
            % Slit correction matrix
            M_deltaslit = [1, Z_slit(ifreq); 0, 1];
            
            % Total transfer matrix
            T_tot = M_deltaslit * (M_s * M_HR * M_s)^M;
            
            % Extract elements
            T11_tot(ifreq) = T_tot(1,1);
            T21_tot(ifreq) = T_tot(2,1);
        end
        
        % Calculate reflection coefficient for this well
        R(:, indn) = (T11_tot - T21_tot * Z0) ./ (T11_tot + T21_tot * Z0);
    end

end

function [L, s_n, R_QRD] = QRDReflectionCoefficient(N,W,design_freq,f_min,f_max,df,rho_0,c_0)

    %% Wavenumber vector
    f_v = (f_min:df:f_max)'; %transposed so that implicit expansion works on line 20
    w=2*pi*f_v;
    k = w./c_0; %wavenumber vector
    
    %% Set parameters
    
    %% QRD depths
    n = (1:N);%vector of desired sequence length
    s_n = mod(n.^2,N);
    lambda_0 = c_0./design_freq; %design wavelength in meters
    L = (s_n.*lambda_0)./(2*N); %well depths
    L(L==0) = 1e-6;
    
    %% Impedance of QWR (no losses)
    z = rho_0.*c_0; %specific acoustic impedance of air at 20degC
    Z_0 = z./W; %characteristic impedance of slit? Careful. This goes to infinity when slit is infinitely thin.
    Zw = 1i.*Z_0.*cot(L.*k);  %array is calculated by implicit expansion
    
    %% Matrices
    R_QRD = ((Zw./Z_0)-1)./((Zw./Z_0)+1);
end
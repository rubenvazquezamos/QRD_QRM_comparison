%Standalone TMM model of N=5 QRD
%
% The main body of the script uses data structures Freq and Geo.
% These are passed to functions and selectively  "unpacked" as needed. 
% e.g. within QRDdepths(), L = Geo.L

clear all
%-----------------------------------------------------------------------%
LatexPreamble() %run LateX preamble
stinson_constants = generateStinsonConstants(); %stinson constants
c_0 = stinson_constants.sound_speed;
rho_0 = stinson_constants.density;

Freq.f_min = 250; % Minimum Freq of interest
Freq.f_max = 4000; % Maximum Freq of interest
Freq.df = 25; % Freq discretization
Freq.dfreq = 500; %design frequency of diffuser [Hz]

Freq.Vector = (Freq.f_min:Freq.df:Freq.f_max);   % Freq vector f_v
Freq.fangular = 2*pi*Freq.Vector; % angular frequency w
Freq.waveno = Freq.fangular./c_0; %wavenumber k
Freq.Nf = numel(Freq.Vector); % Number of frequencies

Geo.N = 5; %prime number used to generate s_n
Geo.W = 6e-2; %width of single well
Geo.Li = 1e-2; %width in between wells
Geo.D = Geo.N*(Geo.W+Geo.Li)+Geo.Li; %implicit
Geo.L = QRDdepths(Geo.N,Freq.dfreq,c_0); %diffuser well depths
Geo.Theta = linspace(0,pi,181);
Geo.numxpts = lcm(Geo.N,100); %number of x points
Geo.x = linspace(0,Geo.D,Geo.numxpts); %x points

%----------------------------------------------------------------------%

%% Reflection coefficient
switches.viscloss = false;
switches.radcorr = true;
R = calculateReflectionCoefficient(Geo,Freq,stinson_constants,switches);
switches.radcorr = false;
R_radcorr = calculateReflectionCoefficient(Geo,Freq,stinson_constants,switches);

%% Diffusion coefficient

input.calcDC = questdlg("calculate diffusion coefficient?");

switch input.calcDC
    case 'Yes'
    % Diffusion coefficient of QRD
    delta.QRD = calculateDiffusionCoefficient(Geo,R,Freq.waveno);
    % Diffusion coefficient with radiation correction
    delta.QRDrc = calculateDiffusionCoefficient(Geo,R_radcorr,Freq.waveno);
    % Flat plane
    Rflat = ones(size(R));
    delta.f = calculateDiffusionCoefficient(Geo,Rflat,Freq.waveno);
    case 'No'
        return
end

input.normDC = questdlg("normalise diffusion coefficient");

switch input.normDC
    case 'Yes'
    % Diffusion coefficient of QRD
    deltan.QRD = (delta.QRD-delta.f)./(1-delta.f);
    % Diffusion coefficient with radiation correction
    deltan.QRDrc = (delta.QRDrc-delta.f)./(1-delta.f);
    case 'No'
end

%% Plotting

figure()
plot(Freq.Vector,delta.QRD,"LineWidth",2)
set(gca,"XScale","log")
hold on
plot(Freq.Vector,delta.QRDrc,"LineWidth",2)
plot(Freq.Vector,delta.f,"LineWidth",2,"LineStyle","--")
ylim([0, 1])
legend(["N=5 QRD","N=5 QRD with rad corr","flat plane"],...
    "Location","southeast")
title("diffusion coefficient (TMM)")
xlabel("Hz")
ylabel("$\delta$")
hold off

if strcmp(input.normDC, 'Yes') %plot normalised diffusion coefficient
    figure()
    plot(Freq.Vector,deltan.QRD,"LineWidth",2)
    hold on
    plot(Freq.Vector,deltan.QRDrc,"LineWidth",2)
    run("comparison_data.m")
    scatter(x_data,y_data)
    ylim([0, 1])
    legend(["N=5 QRD","N=5 QRD with rad corr"],...
        "Location","southeast")
    title("normalised diffusion coefficient (TMM)")
    xlabel("Hz")
    ylabel("$\delta_n$")
    hold off
end

%%---------------------------------------------------------------------------------------------

function [R] = calculateReflectionCoefficient(Geo,Freq,stinson_constants,switches)
%Calculates the reflection coefficient for a single well
%  INPUTS
% 
%   Freq: struct with fields
%       Vector: Frequency points [Hz]      
%       fangular: angular frequency omega
%
%   Geo: struct with fields
%       W: width of wells [m]
%       Li: "fin" width [m]
%       L: vector of well depths [m]
%       N: number of wells and prime number
%
%   stinson_constants: struct with fields
%       c_0: sound speed [m/s]       
%       rho_0: density of free air [kg/m^3]
%       B_0: bulk modulus of free air [Pa]
%       Pr: Prandtl number [dimensionless]
%       eta: kinematic viscosity coefficient
%       gam: specific heat ratio [dimensionless]  
%
%   switches: struct with fields
%       viscloss: [boolean] toggles thin slit losses
%       radcorr: [boolean] toggles radiation correction
%
%  OUTPUTS
%   
%   R: [Nf x N matrix] Complex reflection coefficient
%

    rho_0 = stinson_constants.density;
    c_0 = stinson_constants.sound_speed;
    
    f_v = Freq.Vector;
    w = Freq.fangular;
   
    d = Geo.W+Geo.Li; %total well width
    h = Geo.W; %open well width
    S_0 = Geo.W;
    S_w = Geo.W;

    for indn = 1:Geo.N % Loop over wells
        
        L = Geo.L(indn);      
                
        if switches.viscloss == false 
            c_w = c_0;
        else %incorporate losses into calculation
            c_w = visclossthin(Freq.fangular,Geo.W,stinson_constants);
        end
           
        % Calculate impedances
        Z0 = (rho_0 * c_0) / S_0; % Impedance of free air
        Z_w = (rho_0 .* c_w) ./ S_w; % Impedance of QWR

        % Calculate wavenumber
        k_w = w ./ c_w;

        % Radiation correction
        phi_t = h./d;
        m = 1:1:20; % Summation index
        delta_slit = h * phi_t * sum((sin(m*pi*phi_t).^2) ./ ((m*pi*phi_t).^3));
        
        % characteristic radiation impedance
        Z_deltaQWR = (1i * w * delta_slit * rho_0) / (phi_t * S_0);
        
        if switches.radcorr == false
            Z_deltaQWR = zeros(size(Z_deltaQWR));
        end
        % Waveguide transfer matrix elements
        T11 = cos(k_w * L);
        T12 = 1i * Z_w .* sin(k_w * L);
        T21 = 1i * sin(k_w * L) ./ Z_w;
        T22 = cos(k_w * L);
        
        % Calculate transfer matrices for each frequency
        for ifreq = 1:length(f_v)
            % Waveguide transfer matrix
            M_QWR = [T11(ifreq), T12(ifreq); T21(ifreq), T22(ifreq)];
            
            % Slit correction matrix
            M_deltaQWR = [1, Z_deltaQWR(ifreq); 0, 1];
            
            % Total transfer matrix
            T_tot = M_deltaQWR * M_QWR;
            
            % Extract elements
            T11_tot(ifreq) = T_tot(1,1);
            T21_tot(ifreq) = T_tot(2,1);
        end
        
        % Calculate reflection coefficient for this well
        R(:, indn) = (T11_tot - T21_tot * Z0) ./ (T11_tot + T21_tot * Z0);
    end

end


function delta = calculateDiffusionCoefficient(Geo,R,k)
%Calculates the reflection coefficient for an entire diffuser across
%frequencies
%
% INPUTS
%   Geo: struct with fields
%       N: number of wells and prime number
%       x: vector of points along diffuser surface
%      theta: vector of observer angles
%
%   R: matrix of reflection coefficient values
%   k: vector of wavenumber values
%   
% OUTPUTS
%   delta: vector of diffusion coeffient values
% 
    N = Geo.N;
    theta = Geo.Theta;
    x = Geo.x;


    augind = length(x)./N;
    %augmentation of R matrix so that Fraunhofer Integral can be calculated
    Rbig = repelem(R, 1, augind);    

    % Fraunhofer Integral
    for ifr = 1:length(k)
       Ps(ifr,:) = fftfraunhofer(theta,Rbig(ifr,:),k(ifr),x);
    end
    % CALCULATE DIFFUSION COEFFICIENT
    %-------------------------------------------------------------------------%
    SI = abs(Ps).^2; %SPL is converted to sound intensity.
    n_d = length(theta);
    
    SIsum = sum(SI,2);
    SIsq = sum(SI.^2,2);
    
    delta = (SIsum.^2 - SIsq)./((n_d-1)*(SIsq)); %diffusion coefficient

end


function Ps = fftfraunhofer(theta,Rs,k,x)
    % ======================================
    %
    %   PS = FFTFRAUNHOFER(THETA,RS,K,X)
    %
    %   FFTFRAUNHOFER calculates the far field PS using the Fraunhofer integral
    %   of a surface with reflection RS(X) at wavenumber K and at angles THETA
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

function [c_eff] = visclossthin(w, W,stinson_constants)
%------------------------------------------------------------------------------------
%Calculates the effective acoustic wave speed in a thin slit of 
% width W at angular frequencies w.
%------------------------------------------------------------------------------------
%
% INPUTS
% w angular frequency
% W width of diffuser well
%
%stinson_constants: struct with fields
%   rho_0 density of free air
%   B_0 bulk modulus of free air
%   Pr Prandtl number
%   eta kinematic viscosity
%   gam sepcific heat ratio
%
% OUTPUTS
%c_eff effective wave speed
%
    
    % Unpack required constants
    rho_0 = stinson_constants.density;
    B_0 = stinson_constants.bulk_modulus;
    Pr = stinson_constants.prandtl_no;
    eta = stinson_constants.kinematic_viscosity;
    gam = stinson_constants.specific_heat_ratio;
    
    % Stinson 1991 model
    G_rho = sqrt(-1i*w*rho_0 / eta); %vector
    G_k = sqrt(-1i*w*rho_0*Pr / eta);
    h = W ./ 2; %half widths of wells
    
    % Calculate effective bulk modulus, density and wave speed
    B_eff = B_0 .* ( 1 + (gam-1) .* ( tanh( h .* G_k ) ) ./ ( h .* G_k ) ) .^(-1);
    rho_eff = rho_0 .* ( 1 - ( tanh( h .* G_rho ) ) ./ ( h .* G_rho ) ) .^(-1);
    c_eff = sqrt(B_eff./rho_eff);

end

function L = QRDdepths(N,dfreq,c_0)
% INPUTS
%    N: number of wells and prime number
%    dfreq: design frequency [Hz]
%    c_0: sound speed [m/s]
%
% OUTPUTS
%    L: [1xN double] well depths

    % QRD depths
    n = (0:(N-1));  %vector of desired sequence length
    s_n = mod(n.^2,N); %quadratic residue sequence calculation
    
    lambda_0 = c_0 ./ dfreq; %design wavelength in meters
    L = (s_n.*lambda_0)./(2*N); %well depths
    L(L==0) = 1e-10; % L = 0 not permitted as cot(0)= infty

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


function [stinson_constants] = generateStinsonConstants(~)
    rho_0 = 1.204; %kg/m^3 air density
    B_0 = 1.42e5; %Pa adiabatic bulk modulus
    sp = 1 ; %superficial porosity
    Cv = 0.718 ; %J*K^-1*kg^-1 specific heat stinson_constants volume
    Cp = 1.005; %J*K^-1*kg^-1 specific heat stinson_constants temperature
    mu = 1.825e-5; %N*s/m^2 shear/dynamic viscosity coefficient
    
    %% Derived stinson_constants
    c_0 = sqrt(B_0/rho_0);
    Z_0 = rho_0.*c_0; %Rayls = Pa*s/m impedance of air
    Pr = 0.71;
    eta = mu/rho_0; %kinematic viscosity coeff.
    gam = Cp/Cv; %ratio of specific heats
    
    % Create structure
    stinson_constants.density = rho_0;
    stinson_constants.bulk_modulus = B_0;
    stinson_constants.superficial_porosity = sp;
    stinson_constants.specific_heat_constant_volume = Cv;
    stinson_constants.specific_heat_constant_pressure = Cp;
    stinson_constants.dynamic_viscosity = mu;
    stinson_constants.sound_speed = c_0;
    stinson_constants.impedance = Z_0;
    stinson_constants.prandtl_no = Pr;
    stinson_constants.kinematic_viscosity = eta;
    stinson_constants.specific_heat_ratio = gam;
    
    clear rho_0 B_0 sp Cv Cp mu c_0 Z_0 Pr eta gam

end


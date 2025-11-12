%N=5 QRD diffusion coefficient calculated using the TMM
% Script is meant to be called from within scatteredPressure_main.m

%% Constants
rho = 1.204; %density of air
c = 340; %speed of sound

%% Set parameters
N = 5; %number of wells and prime number
xlabels = 1:N;
fmax = Freq.f_max; % Maximum Freq of interest
fmin = Freq.f_min; % Minimum Freq of interest
df = Freq.df; %frequency step
                                                   
W = Geo.W; %width of well in m.
numxpts = lcm(N,100);
x = linspace(0,Geo.D,numxpts);
theta = linspace(-pi/2,pi/2,181);

% a = Geo.D./2 ; %size of half panel in metres
% W = Geo.W; %width of well in m.
% numxpts = lcm(N,100);
% x = linspace(-a,a,numxpts);
% theta = linspace(-pi/2,pi/2,181);

%% QRD depths
n = (0:(N-1));  %vector of desired sequence length
s_n = mod(n.^2,N);
dfreq = 500; %Hz
lambda_0 = c./dfreq; %design wavelength in meters
Geo.L = (s_n.*lambda_0)./(2*N); %well depths
% Geo.L(Geo.L==0) = 1e-10; % L = 0 not permitted as cot(0)= infty

%% Frequency vectors and QWR geometry
f_v = (fmin:df:fmax)';%frequency vector.
w = f_v.*2.*pi;
k = w./c; %wavenumber vector

%% Impedance of QWR (no losses)
input.radcorr = "Yes";
R = TMM_calc(N,Geo,Freq.f_min, Freq.f_max,Freq.df,rho,c,input);
input.radcorr = "No";
R_radcorr = TMM_calc(N,Geo,Freq.f_min, Freq.f_max,Freq.df,rho,c,input);

%% Function to calculate reflection coefficient
function [R] = TMM_calc(n, Geo, f_min,f_max,df,rho_0,c_0,input)
    f_v = f_min:df:f_max;
    w = 2*pi*f_v;
   
    d = Geo.W+Geo.Li;
    h = Geo.W;
    S_0 = Geo.W;
    S_w = Geo.W;
    
    switch input.radcorr
        case "Yes"
            input.radcorr = 1;
        case "No"
            input.radcorr = 0;
        otherwise
           return;  % e.g. user closed dialog
 end

    for indn = 1:n % Loop over wells
        L = Geo.L(indn);
        % Calculate impedances
        Z0 = (rho_0 * c_0) / S_0; % Impedance of free air
        Z_w = (rho_0 .* c_0) ./ S_w; % Impedance of QWR

        % Calculate wavenumber
        k_w = w ./ c_0;

        % Radiation correction
        phi_t = h./d;
        m = 1:1:20; % Summation index
        delta_slit = h * phi_t * sum((sin(m*pi*phi_t).^2) ./ ((m*pi*phi_t).^3));
        
        % characteristic radiation impedance
        Z_deltaQWR = (1i * w * delta_slit * rho_0) / (phi_t * S_0);
        
        if input.radcorr == 0
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
    

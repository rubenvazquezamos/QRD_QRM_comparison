function R = calculateMetadiffuserReflection(n, geometry, f_min,f_max,df,stinson_constants)
    %calculates the complex reflection coefficient of an acoustic
    %metadiffuser.

    %INPUTS
    %
    % geometry struct with fields
    %   neck.widths/lengths
    %   cavity.widths/lengths
    %   slit.widths/lengths
    %
    %OUTPUTS
    % R Nxnumel(f_v) matrix of complex reflection coefficient values

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



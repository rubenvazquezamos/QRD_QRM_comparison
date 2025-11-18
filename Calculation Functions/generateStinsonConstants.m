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

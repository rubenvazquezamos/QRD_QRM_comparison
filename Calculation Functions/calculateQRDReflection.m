function [L, s_n, R_QRD] = calculateQRDReflection(N,W,design_freq,f_min,f_max,df,rho_0,c_0)

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

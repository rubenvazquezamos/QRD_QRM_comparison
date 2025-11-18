function delta = calculateDiffusionCoefficient(Ps,theta,ri)
    % INPUTS:
    %   Ps: [NxM] matrix of scattered pressure across N wavenumbers 
    %       and N observer angles
    %   theta: vector of observer angles
    %   ri: scalar restriction on extreme observer angles
    %      eg: ri = 2 gives angle range -88 to 88 deg.
    % OUTPUT:
    %   delta: vector of diffusion coefficient values
    %-------------------------------------------------------------------------%
   
    if ri ~= 0
      theta = Geo.Theta(ri:end-ri); %restrict extreme observer angles
      Ps = Ps(:,ri:end-ri);
    end

    SI = abs(Ps).^2; %SPL is converted to sound intensity.
    n_d = length(theta);
    
    SIsum = sum(SI,2);
    SIsq = sum(SI.^2,2);
    
    delta = (SIsum.^2 - SIsq)./((n_d-1)*(SIsq)); %diffusion coefficient

end


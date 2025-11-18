function delta_n = normaliseDiffusionCoefficient(delta, deltaf)
    %INPUTS
    % delta: diffusion coefficient of reflector
    % deltaf: diffusion coeffficient of flat plane
    %
    %OUTPUT
    % delta_n: normalised diffusion coefficient

    delta_n = (delta - deltaf)./(1-deltaf);
end
function Rbig = discretizeDiffuserReflection(N,x,R)
    % INPUTS
    %   N: number of wells and prime number
    %   x: vector of points along diffuser surface.
    %   R: matrix of reflection coefficient values.
    % 
    % OUTPUTS
    %   Rbig: matrix of spatially discretized reflection coefficient values.
            
        augmentation_index = length(x)./N;
        Rbig = repelem(R, 1, augmentation_index);    
end

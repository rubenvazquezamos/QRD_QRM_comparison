function Ps = fftfraunhofer(theta,Rs,k,x)
    % ======================================
    %
    %   PS = FFTFRAUNHOFER(THETA,RS,K,X)
    %
    %   FFTFRAUNHOFER calculates the far field PS using the Fraunhofer integral
    %   of a surface with reflection RS(X) at wavenumber k and at angles
    %   theta
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


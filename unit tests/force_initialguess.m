% unit test to check conversion of geometry structure into geometry matrix
% useful when forcing an initial guess

initialguess = cell2mat(flattenStruct2Cell(load('optgeo1.mat')));
flatinitialguess = initialguess';
flatinitialguess = flatinitialguess(:); %ballpark initial guess

checkinitialguess = unpackgeometry(flatinitialguess,5);
geomchecker(checkinitialguess,5,12e-3); %check it's been flattened correctly

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
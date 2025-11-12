%% Optimisation setup
% Inequality constraints
% [w_n,l_n,w_c,l_c,h,a_y]
A = [
    0 0  1 0 0 -1;
    1 0 -1 0 0  0;
                  ];
b = [0; 0];

% Equality constraints
% [w_n,l_n,w_c,l_c,h,a_y]
Aeq = [
        0 0 0 0 0 1;
        0 0 1 0 0 0
        0 1  0 1 1  0
                    ];

beq = [L; L-e, ; a_x];

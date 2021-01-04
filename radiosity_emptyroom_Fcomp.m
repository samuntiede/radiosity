% Calculate the F matrix for the cube room with all six square walls.
% See the video https://youtu.be/krIVZvzlxUQ
%
% Samuli Siltanen January 2021

%% Preliminaries

% Choose mosaic resolution. Each wall, floor anf roof surface will be
% divided into n x n square-shaped patches. See 
n = 10;

% Choose integration quadrature parameter. Integrals over pixels are
% implemented as midpoint rule using qn x qn sub-pixels of constant size
qn = 3;

% Construct centerpoints of the mosaic tiles.
% The letter d denotes the length of the side of a pixel.
d        = 2/n;
tmp      = -1-d/2 + [1:n]*d;

% Initialize centerpoint coordinate matrices
Xmat = zeros(n^2,6);
Ymat = zeros(n^2,6);
Zmat = zeros(n^2,6);

% Construct the centerpoints for all the tiles in all the six walls.
% The ordering of the five walls below fixes the indexing of all the tiles
% using just one number running from 1 to 6*(n^2).

% The back wall (1)
[X,Z]     = meshgrid(tmp);
Xmat(:,1) = X(:);
Zmat(:,1) = Z(:);
Ymat(:,1) = ones(n^2,1);

% Roof (2)
[X,Y] = meshgrid(tmp);
Xmat(:,2) = X(:);
Ymat(:,2) = Y(:);
Zmat(:,2) = ones(n^2,1);

% Floor (3)
Xmat(:,3) = X(:);
Ymat(:,3) = Y(:);
Zmat(:,3) = -ones(n^2,1);

% Right-hand-side wall (4)
[Y,Z] = meshgrid(tmp);
Ymat(:,4) = Y(:);
Zmat(:,4) = Z(:);
Xmat(:,4) = ones(n^2,1);

% Left-hand-side wall (5)
Ymat(:,5) = Y(:);
Zmat(:,5) = Z(:);
Xmat(:,5) = -ones(n^2,1);

% Front wall (6), invisible in renderings
[X,Z]     = meshgrid(tmp);
Xmat(:,6) = X(:);
Zmat(:,6) = Z(:);
Ymat(:,6) = -ones(n^2,1);

% Formula for view factor between square-shaped pixels sharing an edge.
% From Cohen&Wallace: Radiosity and realistic image synthesis
% (Academic Press Professional 1993), Figure 4.4
shared_edge_F = (2*atan(1)-sqrt(2)*atan(1/sqrt(2))+.25*log(3/4))/pi;

% Quadrature points and weights for integrating over a square of size d x d
% centered at the origin
tt = [1:qn]/qn*d - .5*d/qn - d/2;
[q1,q2] = meshgrid(tt);
qw = (d/qn)^4; % Area of quadrature pixel, squared, serves as the weight


%% Form the geometrical view factor matrix F.
% See http://en.wikipedia.org/wiki/View_factor for details of computation.

% Initialize the matrix
F = zeros(6*n^2);

% From the roof (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the roof
        pjjj = [Xmat(jjj,2);Ymat(jjj,2);Zmat(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(iii,n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the roof pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors roof->back done (1/15)')

% From the floor (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the floor
        pjjj = [Xmat(jjj,3);Ymat(jjj,3);Zmat(jjj,3)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(iii,2*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the floor pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(3));
                    cosiii = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,2*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors floor->back done (2/15)')

% From the right-hand-side wall (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the roof
        pjjj = [Xmat(jjj,4);Ymat(jjj,4);Zmat(jjj,4)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(iii,3*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the right wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(1));
                    cosiii = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,3*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors right->back done (3/15)')

% From the left-hand-side wall (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the roof
        pjjj = [Xmat(jjj,5);Ymat(jjj,5);Zmat(jjj,5)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(iii,4*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the left wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(1));
                    cosiii = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,4*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors left->back done (4/15)')


% From the front wall (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the front wall
        pjjj = [Xmat(jjj,6);Ymat(jjj,6);Zmat(jjj,6)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(iii,5*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the front wall pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,5*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors front->back done (5/15)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the floor (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the floor
        pjjj = [Xmat(jjj,3);Ymat(jjj,3);Zmat(jjj,3)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(n^2+iii,2*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the floor pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(3));
                    cosiii = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(n^2+iii,2*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors floor->roof done (6/15)')

% From the right-hand wall (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the right wall
        pjjj = [Xmat(jjj,4);Ymat(jjj,4);Zmat(jjj,4)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(n^2+iii,3*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the right wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(1));
                    cosiii = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(n^2+iii,3*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors right->roof done (7/15)')

% From the left-hand-side wall (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the left wall
        pjjj = [Xmat(jjj,5);Ymat(jjj,5);Zmat(jjj,5)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(n^2+iii,4*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the left wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(1));
                    cosiii = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(n^2+iii,4*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors left->roof done (8/15)')


% From the front wall (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the front wall
        pjjj = [Xmat(jjj,6);Ymat(jjj,6);Zmat(jjj,6)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(n^2+iii,5*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the front wall pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(3));
                    cosjjj = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(n^2+iii,5*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors front->roof done (9/15)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the right-hand-side wall (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the right-hand-side wall
        pjjj = [Xmat(jjj,4);Ymat(jjj,4);Zmat(jjj,4)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(2*n^2+iii,3*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the right wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosjjj = abs(tmp2(1));
                    cosiii = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(2*n^2+iii,3*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors right->floor done (10/15)')

% From the left-hand-side wall (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the left wall
        pjjj = [Xmat(jjj,5);Ymat(jjj,5);Zmat(jjj,5)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(2*n^2+iii,4*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the left wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(3));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(2*n^2+iii,4*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors left->floor done (11/15)')

% From the front wall (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the front wall
        pjjj = [Xmat(jjj,6);Ymat(jjj,6);Zmat(jjj,6)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(2*n^2+iii,5*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the front wall pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(3));
                    cosjjj = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(2*n^2+iii,5*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors front->floor done (12/15)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the left-hand-side wall (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the left wall
        pjjj = [Xmat(jjj,5);Ymat(jjj,5);Zmat(jjj,5)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(3*n^2+iii,4*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the left wall pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(1));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(3*n^2+iii,4*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors left->right done (13/15)')

% From the front wall (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the front wall
        pjjj = [Xmat(jjj,6);Ymat(jjj,6);Zmat(jjj,6)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(3*n^2+iii,5*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the front wall pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(1));
                    cosjjj = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(3*n^2+iii,5*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors front->right done (14/15)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the front wall (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:n^2
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in the front wall
        pjjj = [Xmat(jjj,6);Ymat(jjj,6);Zmat(jjj,6)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(4*n^2+iii,5*n^2+jjj) = shared_edge_F;
        else % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the left wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the front wall pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(1));
                    cosjjj = abs(tmp2(2));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(4*n^2+iii,5*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors front->right done (15/15)')


% Use symmetry to finish the construction of F. F is symmetric
% since all the pixels in our model have equal size.
F = F+F.';

% Check the matrix F. The row sums should all be one
figure(10)
clf
plot(sum(F))
title('Check plot: all values should ideally be one')


% Save matrix to disc
save data/F_emptyroom F n qn d Xmat Ymat Zmat 



% Further testing for building a radiosity lighting method for virtual
% spaces. Here we implement a table
%
% Samuli Siltanen January 2021

%% Preliminaries

% Choose mosaic resolution; n must be even
halfn = 12;
n = 2*halfn;
n_wall = 6*n^2;
n_table = 2*halfn^2+8*halfn;

% Define the location of the table
t_left  = -.5; % Must be strictly between -1 and 0
t_front = -.5; % Must be strictly between -1 and 0
t_top   = -.70; % Must be strictly between -1+2*d and 1

% Choose integration quadrature parameter. Integrals over pixels are
% implemented as midpoint rule using qn x qn sub-pixels of constant size
qn = 3;

% Construct centerpoints of the mosaic tiles.
% The letter d denotes the length of the side of a pixel.
d        = 2/n;
tmp      = -1-d/2 + [1:n]*d;

% Formula for view factor between square-shaped pixels sharing an edge.
% From Cohen&Wallace: Radiosity and realistic image synthesis
% (Academic Press Professional 1993), Figure 4.4
shared_edge_F = (2*atan(1)-sqrt(2)*atan(1/sqrt(2))+.25*log(3/4))/pi;

% Quadrature points and weights for integrating over a square of size d x d
% centered at the origin
tt = [1:qn]/qn*d - .5*d/qn - d/2;
[q1,q2] = meshgrid(tt);
qw = (d/qn)^4; % Area of quadrature pixel, squared, serves as the weight


%% Construct the centerpoints for all the tiles in all the six walls.

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


%% Construct the centerpoints for all the tiles in the table
% The table is chosen to be two pixels thick

% Initialize centerpoint coordinate matrices for the table top and bottom
Xmat2 = zeros(halfn^2,2);
Ymat2 = zeros(halfn^2,2);
Zmat2 = zeros(halfn^2,2);
tmp1 = -d/2 + [1:halfn]*d;
% Table top surface (7)
[X,Y] = meshgrid(tmp1,tmp1);
Xmat2(:,1) = t_left+X(:);
Ymat2(:,1) = t_front+Y(:);
Zmat2(:,1) = t_top*ones(halfn^2,1);
% Table bottom surface (8)
Xmat2(:,2) = t_left+X(:);
Ymat2(:,2) = t_front+Y(:);
Zmat2(:,2) = (t_top-2*d)*ones(halfn^2,1);

% Initialize centerpoint coordinate matrices for the table sides
Xmat3 = zeros(2*halfn,4);
Ymat3 = zeros(2*halfn,4);
Zmat3 = zeros(2*halfn,4);
% Back side of the table (9)
Xmat3(:,1) = t_left+[tmp1(:);tmp1(:)];
Ymat3(:,1) = t_front+halfn*d*ones(2*halfn,1);
Zmat3(:,1) = t_top-[d/2*ones(halfn,1);3*d/2*ones(halfn,1);];
% Right side of the table (10)
Xmat3(:,2) = t_left+halfn*d*ones(2*halfn,1);
Ymat3(:,2) = t_front+[tmp1(:);tmp1(:)];
Zmat3(:,2) = Zmat3(:,1);
% Left side of the table (11)
Xmat3(:,3) = t_left*ones(2*halfn,1);
Ymat3(:,3) = Ymat3(:,2);
Zmat3(:,3) = Zmat3(:,2);
% Front side of the table (12)
Xmat3(:,4) = Xmat3(:,1);
Ymat3(:,4) = t_front*ones(2*halfn,1);
Zmat3(:,4) = Zmat3(:,1);



%% Form the geometrical view factor matrix F.
% See http://en.wikipedia.org/wiki/View_factor for details of computation.

% Initialize the matrix
F = zeros(n_wall+n_table);

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
        tmp2    = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors roof->back wall done (1/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors floor->back wall done (2/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's right edge?
        tmp_right = pjjj+(t_left+halfn*d-pjjj(1))*tmp2/tmp2(1);
        if tmp_right(3)>t_top-2*d & tmp_right(3)<t_top & tmp_right(2)>t_front & tmp_right(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors right->back wall done (3/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's left edge?
        tmp_left = pjjj+(t_left-pjjj(1))*tmp2/tmp2(1);
        if tmp_left(3)>t_top-2*d & tmp_left(3)<t_top & tmp_left(2)>t_front & tmp_left(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors left->back wall done (4/66)')

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
        tmp2    = difvec0/r0;
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's front edge?
        tmp_front = pjjj+(t_front-pjjj(2))*tmp2/tmp2(2);
        if tmp_front(3)>t_top-2*d & tmp_front(3)<t_top & tmp_front(1)>t_left & tmp_front(1)<(t_left+halfn*d)
            visible = 0;
        end
        % Does the connecting line cross the table's back edge?
        tmp_back = pjjj+(t_front+halfn*d-pjjj(2))*tmp2/tmp2(2);
        if tmp_back(3)>t_top-2*d & tmp_back(3)<t_top & tmp_back(1)>t_left & tmp_back(1)<(t_left+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
            % Edge never shared: integrate for F using quadrature
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
disp('Geometric view factors front->back wall done (5/66)')

% From the table top (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in table top
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table top pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
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
            F(iii,6*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tabletop->back wall done (6/66)')

% From the table bottom (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the table bottom
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(3)<pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table bottom pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
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
            F(iii,6*n^2+halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tablebottom->back wall done (7/66)')

% From the table back side (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the table back side
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        % Edge never shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the back wall pixel
                qpiii = piii;
                qpiii(1) = qpiii(1)+q1(kkk);
                qpiii(3) = qpiii(3)+q2(kkk);
                % Quadrature point in the table back side pixel
                qpjjj = pjjj;
                qpjjj(1) = qpjjj(1)+q1(lll);
                qpjjj(2) = qpjjj(2)+q2(lll);
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
        F(iii,6*n^2+2*halfn^2+jjj) = viewfactor;
    end
end
disp('Geometric view factors tableback->back wall done (8/66)')

% From the table right side (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in table right side
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)>pjjj(1) % visibility check
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table right side pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tableright->back wall done (9/66)')

% From the table left side (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in table left side
        pjjj = [Xmat3(jjj,3);Ymat3(jjj,3);Zmat3(jjj,3)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)<pjjj(1) % visibility check
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the back wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table left side pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(iii,6*n^2+2*halfn^2+4*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tableleft->back wall done (10/66)')

% From the table front surface (jjj) to the back wall (iii)
% There is no visibility at all
disp('Geometric view factors tableleft->back wall done (11/66)')







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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
            % Edge never shared: integrate for F using quadrature
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
disp('Geometric view factors floor->roof done (12/66)')

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
        tmp2    = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors right->roof done (13/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors left->roof done (14/66)')

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
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors front->roof done (15/66)')

% From the table top (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        % Check if the two pixels share an edge
        if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
            % Calculate element of F analytically
            F(n^2+iii,6*n^2+jjj) = shared_edge_F;
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
                    % Quadrature point in the cubicle-right pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(3));
                    cosjjj = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(n^2+iii,6*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tabletop->roof done (16/66)')

% From the table bottom (jjj) to the roof (iii)
% No visibility
disp('Geometric view factors tablebottom->roof done (17/66)')

% From the table back side (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the table back edge
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        % Edge not ever shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the roof pixel
                qpiii = piii;
                qpiii(1) = qpiii(1)+q1(kkk);
                qpiii(2) = qpiii(2)+q2(kkk);
                % Quadrature point in the cubicle-top pixel
                qpjjj = pjjj;
                qpjjj(1) = qpjjj(1)+q1(lll);
                qpjjj(2) = qpjjj(2)+q2(lll);
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
        F(n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
    end
end
disp('Geometric view factors tableback->roof done (18/66)')

% From the table right side (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)>pjjj(1) % Visibility check
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the table-right pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
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
            F(n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tableright->roof done (19/66)')

% From the table left side (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,3);Ymat3(jjj,3);Zmat3(jjj,3)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)>pjjj(1) % Visibility check
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the cubicle-top pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
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
            F(n^2+iii,6*n^2+2*halfn^2+4*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tableleft->roof done (20/66)')

% From the table front side (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the table front
        pjjj = [Xmat3(jjj,4);Ymat3(jjj,4);Zmat3(jjj,4)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(2)<pjjj(2) % Visibility check
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the table front pixel
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
            F(n^2+iii,6*n^2+2*halfn^2+6*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tablefront->roof done (21/66)')






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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors right wall->floor done (22/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors left->floor done (23/66)')

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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        tmp3 = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
        if tmp3(1)>t_left & tmp3(1)<(t_left+halfn*d) & tmp3(2)>t_front & tmp3(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        tmp4 = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
        if tmp4(1)>t_left & tmp4(1)<(t_left+halfn*d) & tmp4(2)>t_front & tmp4(2)<(t_front+halfn*d) % Does the connecting line cross the table?
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors front->floor done (24/66)')

% From the tabletop (jjj) to the floor (iii)
% No visibility
disp('Geometric view factors tabletop->floor done (25/66)')

% From the table bottom (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the table bottom
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        % Edge never shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the floor pixel
                qpiii = piii;
                qpiii(1) = qpiii(1)+q1(kkk);
                qpiii(2) = qpiii(2)+q2(kkk);
                % Quadrature point in the table bottom pixel
                qpjjj = pjjj;
                qpjjj(2) = qpjjj(2)+q1(lll);
                qpjjj(3) = qpjjj(3)+q2(lll);
                % Vector connecting the quadrature points
                difvec = qpiii-qpjjj;
                r      = norm(difvec);
                tmp2   = difvec/r; % Unit direction vector
                cosiii = abs(tmp2(3));
                cosjjj = abs(tmp2(3));
                % Evaluate integrand
                intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
            end
        end
        % Calculate element of F
        viewfactor = qw*sum(sum(intgndmat))/d^2;
        F(2*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
    end
end
disp('Geometric view factors tablebottom->floor done (26/66)')

% From the table back edge (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % From the table back (jjj) to the floor (iii)
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the table back edge
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(2)>pjjj(2)
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the roof pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the table back edge pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
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
            F(2*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tableback->floor (27/66)')

% From the table right edge (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in table right edge
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)>pjjj(1) % visibility check
            % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the  table right edge pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
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
            F(2*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table right edge->floor done (28/66)')

% From the table left edge (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,3);Ymat3(jjj,3);Zmat3(jjj,3)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(1)<pjjj(1) % visibility check
            % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the  table left edge pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
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
            F(2*n^2+iii,6*n^2+2*halfn^2+4*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table left edge->floor done (29/66)')

% From the table front edge (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the table front
        pjjj = [Xmat3(jjj,4);Ymat3(jjj,4);Zmat3(jjj,4)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(2)<pjjj(2) % visibility check
            % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the floor pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(2) = qpiii(2)+q2(kkk);
                    % Quadrature point in the  table front edge pixel
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
            F(2*n^2+iii,6*n^2+2*halfn^2+6*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table front edge->floor done (30/66)')






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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's left edge?
        tmp_left = pjjj+(t_left-pjjj(1))*tmp2/tmp2(1);
        if tmp_left(3)>t_top-2*d & tmp_left(3)<t_top & tmp_left(2)>t_front & tmp_left(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Does the connecting line cross the table's right edge?
        tmp_right = pjjj+(t_left+halfn*d-pjjj(1))*tmp2/tmp2(1);
        if tmp_right(3)>t_top-2*d & tmp_right(3)<t_top & tmp_right(2)>t_front & tmp_right(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
            % Edge not shared: integrate for F using quadrature
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
disp('Geometric view factors left->right-hand-side wall done (31/66)')

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
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's right edge?
        tmp_right = pjjj+(t_left+halfn*d-pjjj(1))*tmp2/tmp2(1);
        if tmp_right(3)>t_top-2*d & tmp_right(3)<t_top & tmp_right(2)>t_front & tmp_right(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors front->right-hand-side wall done (32/66)')

% From the table top (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in table top
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the tabletop pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(1));
                    cosjjj = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(3*n^2+iii,6*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tabletop->right-hand-side wall done (33/66)')

% From the table bottom (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the table bottom
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(3)<pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge not shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table bottom pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(1));
                    cosjjj = abs(tmp2(3));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(3*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors tablebottom->right-hand-side wall done (34/66)')

% From the table back (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in table back
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(2)>pjjj(2) % visibility check
            % Edge not ever shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table back pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
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
            F(3*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table back->right-hand-side wall (35/66)')

% From the table right (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the table right
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Edge never shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the right-hand-side wall pixel
                qpiii = piii;
                qpiii(2) = qpiii(2)+q1(kkk);
                qpiii(3) = qpiii(3)+q2(kkk);
                % Quadrature point in the table right pixel
                qpjjj = pjjj;
                qpjjj(1) = qpjjj(1)+q1(lll);
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
        F(3*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
    end
end
disp('Geometric view factors table right->right-hand-side wall done (36/66)')

% From the table left (jjj) to the right-hand-side wall (iii)
% No visibility
disp('Geometric view factors table left->right-hand-side wall done (37/66)')

% From the table front (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the table front
        pjjj = [Xmat3(jjj,4);Ymat3(jjj,4);Zmat3(jjj,4)];
        if pjjj(2)>piii(2)
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the right-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table front pixel
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
            F(3*n^2+iii,6*n^2+2*halfn^2+6*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table front->right-hand-side wall done (38/66)')







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
        tmp2   = difvec0/r0;
        % Check visibility
        visible = 1;
        if abs(piii(3)-pjjj(3))>d/2 % Line connecting the points NOT horizontal
            % Does the connecting line cross the table top surface?
            tmp_top = pjjj+(t_top-pjjj(3))*tmp2/tmp2(3);
            if tmp_top(1)>t_left & tmp_top(1)<(t_left+halfn*d) & tmp_top(2)>t_front & tmp_top(2)<(t_front+halfn*d) 
                visible = 0;
            end
            % Does the connecting line cross the table bottom surface?
            tmp_bottom = pjjj+(t_top-2*d-pjjj(3))*tmp2/tmp2(3);
            if tmp_bottom(1)>t_left & tmp_bottom(1)<(t_left+halfn*d) & tmp_bottom(2)>t_front & tmp_bottom(2)<(t_front+halfn*d) 
                visible = 0;
            end
        end
        % Does the connecting line cross the table's left edge?
        tmp_left = pjjj+(t_left-pjjj(1))*tmp2/tmp2(1);
        if tmp_left(3)>t_top-2*d & tmp_left(3)<t_top & tmp_left(2)>t_front & tmp_left(2)<(t_front+halfn*d)
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors front->left-hand-side wall done (39/66)')

% From the table top (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in table top
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(4*n^2+iii,6*n^2+halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the table top pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
                        qpjjj(3) = qpjjj(3)+q2(lll);
                        % Vector connecting the quadrature points
                        difvec = qpiii-qpjjj;
                        r      = norm(difvec);
                        tmp2   = difvec/r; % Unit direction vector
                        cosiii = abs(tmp2(1));
                        cosjjj = abs(tmp2(3));
                        % Evaluate integrand
                        intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                    end
                end
                % Calculate element of F
                viewfactor = qw*sum(sum(intgndmat))/d^2;
                F(4*n^2+iii,6*n^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors tabletop->left-hand-side wall done (40/66)')

% From the table bottom (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in table bottom
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(3)<pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(4*n^2+iii,6*n^2+halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the table bottom pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
                        qpjjj(3) = qpjjj(3)+q2(lll);
                        % Vector connecting the quadrature points
                        difvec = qpiii-qpjjj;
                        r      = norm(difvec);
                        tmp2   = difvec/r; % Unit direction vector
                        cosiii = abs(tmp2(1));
                        cosjjj = abs(tmp2(3));
                        % Evaluate integrand
                        intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                    end
                end
                % Calculate element of F
                viewfactor = qw*sum(sum(intgndmat))/d^2;
                F(4*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors table bottom->left-hand-side wall done (41/66)')

% From the table back (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in the table back
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(2)>pjjj(2) % visibility check
            % Edge not ever shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the left-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table back pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
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
            F(4*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table back->left-hand-side wall (42/66)')

% From the table right (jjj) to the left-hand-side wall (iii)
% No visibility
disp('Geometric view factors table right ->left-hand-side wall done (43/66)')

% From the table left (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in table left
        pjjj = [Xmat3(jjj,3);Ymat3(jjj,3);Zmat3(jjj,3)];
        % Edge never shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the left-hand-side wall pixel
                qpiii = piii;
                qpiii(2) = qpiii(2)+q1(kkk);
                qpiii(3) = qpiii(3)+q2(kkk);
                % Quadrature point in the table left pixel
                qpjjj = pjjj;
                qpjjj(1) = qpjjj(1)+q1(lll);
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
        F(4*n^2+iii,6*n^2+2*halfn^2+4*halfn+jjj) = viewfactor;
    end
end
disp('Geometric view factors table left->left-hand-side wall done (44/66)')


% From the table front (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in table front
        pjjj = [Xmat3(jjj,4);Ymat3(jjj,4);Zmat3(jjj,4)];
        if piii(2)<pjjj(2)
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the left-hand-side wall pixel
                    qpiii = piii;
                    qpiii(2) = qpiii(2)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table front pixel
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
            F(4*n^2+iii,6*n^2+2*halfn^2+6*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table front->left-hand-side wall done (45/66)')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From the table top (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in table top
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the front wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table top pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
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
            F(5*n^2+iii,6*n^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table top->front wall done (46/66)')

% From the table bottom (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in table bottom
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(3)<pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the front wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table bottom pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
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
            F(5*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table bottom->front wall done (47/66)')

% From table back (jjj) to the front wall (iii)
% No visibility
disp('Geometric view factors table back->front wall done (48/66)')

% From table right (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in the table right edge
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        if piii(1)>pjjj(1) % visibility check
            % Edge not ever shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the front wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table right pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(5*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table right->front wall (49/66)')

% From the table left (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in table left
        pjjj = [Xmat3(jjj,3);Ymat3(jjj,3);Zmat3(jjj,3)];
        if piii(1)<pjjj(1)
            % Edge never shared: integrate for F using quadrature
            % Initialize matrix of integrand values at quadrature points
            intgndmat = zeros(qn^2,qn^2);
            % Double loop over four-dimensional quadrature
            for kkk = 1:qn^2
                for lll = 1:qn^2
                    % Quadrature point in the front wall pixel
                    qpiii = piii;
                    qpiii(1) = qpiii(1)+q1(kkk);
                    qpiii(3) = qpiii(3)+q2(kkk);
                    % Quadrature point in the table left pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(3) = qpjjj(3)+q2(lll);
                    % Vector connecting the quadrature points
                    difvec = qpiii-qpjjj;
                    r      = norm(difvec);
                    tmp2   = difvec/r; % Unit direction vector
                    cosiii = abs(tmp2(2));
                    cosjjj = abs(tmp2(1));
                    % Evaluate integrand
                    intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
                end
            end
            % Calculate element of F
            viewfactor = qw*sum(sum(intgndmat))/d^2;
            F(5*n^2+iii,6*n^2+2*halfn^2+4*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors table left->front wall done (50/66)')

% From the table front (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in table front
        pjjj = [Xmat3(jjj,4);Ymat3(jjj,4);Zmat3(jjj,4)];
        % Edge never shared: integrate for F using quadrature
        % Initialize matrix of integrand values at quadrature points
        intgndmat = zeros(qn^2,qn^2);
        % Double loop over four-dimensional quadrature
        for kkk = 1:qn^2
            for lll = 1:qn^2
                % Quadrature point in the front wall pixel
                qpiii = piii;
                qpiii(1) = qpiii(1)+q1(kkk);
                qpiii(3) = qpiii(3)+q2(kkk);
                % Quadrature point in the table front pixel
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
        F(5*n^2+iii,6*n^2+2*halfn^2+6*halfn+jjj) = viewfactor;
    end
end
disp('Geometric view factors table front->front wall done (51/66)')

% The sides of the table have no visibility to each other
disp('The rest of geometric view factors are zero')


% Use symmetry to finish the construction of F. F is symmetric
% since all the pixels in our model have equal size.
F = F+F.';

% Check the matrix F. The row sums should all be one
figure(10)
clf
plot(sum(F))
% figure(11)
% clf
% spy(F)


% Save matrix to disc
save data/F_table F n qn d Xmat Ymat Zmat Xmat2 Ymat2 Zmat2 Xmat3 Ymat3 Zmat3 halfn n_wall n_table t_top t_front t_left
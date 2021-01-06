% Calculate the F matrix for the cube room with all six square walls and an
% internal dividing wall. 
% See the video https://youtu.be/krIVZvzlxUQ
%
% Samuli Siltanen January 2021

%% Preliminaries

% Choose mosaic resolution; n must be even. That's why we actually set 
% choose halfn and set n = 2*halfn;
%
% Each wall, floor anf roof surface will be divided into n x n square-shaped 
% patches. Furthermore, the thickness of the dividing wall is two patches
% and its sides have halfn x halfn patches. 
%
% Value halfn=10 is good for making quick tests. 
% Bigger n gives more details, but also the calculation and
% rendering becomes more heavy for the computer. For my computational
% equipment halfn=24 is about the biggest I can use without running out of
% memory. 
halfn = 10;
n = 2*halfn;
n_wall = 6*n^2;
n_cubicle = 2*halfn^2+4*halfn;

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


%% Construct the centerpoints for all the tiles in the "cubicle wall"
% The cubicle wall is chosen to be two pixels thick

% Small gap between cubicle and other surfaces, for avoiding geometric drawing errors
smgap = .05*d;

% Initialize centerpoint coordinate matrices for the two walls
Xmat2 = zeros(halfn^2,2);
Ymat2 = zeros(halfn^2,2);
Zmat2 = zeros(halfn^2,2);
tmp1 = -d/2 + [1:halfn]*d;
% Right-hand-side cubicle wall (7)
[Y,Z] = meshgrid(tmp1,tmp1-1);
Ymat2(:,1) = Y(:)-smgap;
Zmat2(:,1) = Z(:)+smgap;
Xmat2(:,1) = (-1/2+2*d)*ones(halfn^2,1);
% Left-hand-side cubicle wall (8)
Ymat2(:,2) = Y(:)-smgap;
Zmat2(:,2) = Z(:)+smgap;
Xmat2(:,2) = (-1/2)*ones(halfn^2,1);


% Initialize centerpoint coordinate matrices for the top and front
Xmat3 = zeros(2*halfn,2);
Ymat3 = zeros(2*halfn,2);
Zmat3 = zeros(2*halfn,2);
% Top surface of the cubicle wall (9)
Xmat3(:,1) = [(-1/2+d/2)*ones(halfn,1);(-1/2+3*d/2)*ones(halfn,1);];
Ymat3(:,1) = [tmp1(:);tmp1(:)]-smgap;
Zmat3(:,1) = zeros(2*halfn,1)+smgap;
% Front surface of the cubicle wall (10)
Xmat3(:,2) = [(-1/2+d/2)*ones(halfn,1);(-1/2+3*d/2)*ones(halfn,1);];
Ymat3(:,2) = zeros(2*halfn,1)-smgap;
Zmat3(:,2) = [tmp1(:)-1;tmp1(:)-1]+smgap;



%% Form the geometrical view factor matrix F.
% See http://en.wikipedia.org/wiki/View_factor for details of computation.

% Initialize the matrix
F = zeros(n_wall+n_cubicle);

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
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0 % Back wall pixels behind cubicle
            visible = 0;
        end
        % If the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
        end
        % If the back wall centerpoint is within cubicle wall, there is no visibility
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0
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
disp('Geometric view factors roof->back done (1/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0 % Back wall pixels behind cubicle
            visible = 0;
        end
        if abs(pjjj(1)-(-1/2+d))<d & pjjj(2)>0 % Floor pixels behind cubicle
            visible = 0;
        end
        % If the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
        end
        % If the two centerpoints both have x-coordinates within cubicle wall,
        % there is invisibility
        if abs(piii(1)-(-1/2+d))<d & abs(pjjj(1)-(-1/2+d))<d
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
disp('Geometric view factors floor->back done (2/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0 % Back wall pixels behind cubicle
            visible = 0;
        end
        % Only if the back wall patch is on the left of the cubicle wall,
        % there is a chance of invisibility
        if piii(1)<(-1/2)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors right->back done (3/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0 % Back wall pixels behind cubicle
            visible = 0;
        end
        % Only if the back wall patch is on the right of the cubicle wall,
        % there is a chance of invisibility
        if piii(1)>(-1/2+2*d)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors left->back done (4/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(3)<0 % Back wall pixels behind cubicle
            visible = 0;
        end
        % Only if the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors front->back done (5/45)')

% From the right-hand-side cubicle wall (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(1)>pjjj(1) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(iii,6*n^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-right pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
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
                F(iii,6*n^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleright->back done (6/45)')

% From the left-hand-side cubicle wall (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the left-hand-side cubicle wall
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(1)<pjjj(1) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(iii,6*n^2+halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-left pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
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
                F(iii,6*n^2+halfn^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleleft->back done (7/45)')

% From the cubicle top surface (jjj) to the back wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the back wall
        piii = [Xmat(iii,1);Ymat(iii,1);Zmat(iii,1)];
        % Centerpoint of the current pixel in the cubicle top
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(iii,6*n^2+2*halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-top pixel
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
                F(iii,6*n^2+2*halfn^2+jjj) = viewfactor;
                
            end
        end
    end
end
disp('Geometric view factors cubicletop->back done (8/45)')
    
% From the cubicle front surface (jjj) to the back wall (iii)
% There is no visibility at all
disp('Geometric view factors cubiclefront->back done (9/45)')







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
        if abs(pjjj(1)-(-1/2+d))<d & pjjj(2)>0 % Floor pixels behind cubicle
            visible = 0;
        end
        % If the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
        end
                % If the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
        end
        % If the two centerpoints both have x-coordinates within cubicle wall,
        % there is possible invisibility
        if abs(piii(1)-(-1/2+d))<d & abs(pjjj(1)-(-1/2+d))<d & piii(2)>0 & pjjj(2)> -piii(2)
            visible = 0;
        end
        if abs(pjjj(1)-(-1/2+d))<d & pjjj(2)>0
            visible = 0;
        end
        % Calculate nonzero element of F only if visibility holds
        if visible
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
end
disp('Geometric view factors floor->roof done (10/45)')

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
        % Check visibility: no visibility check needed as the shadow effect
        % is very weak between roof and back wall
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
disp('Geometric view factors right->roof done (11/45)')

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
        % Only if the roof patch is on the right of the cubicle wall,
        % there is a chance of invisibility
        if piii(1)>(-1/2+2*d)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors left->roof done (12/45)')


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
        % Check visibility: no visibility check needed as the shadow effect
        % is very weak between roof and back wall
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
disp('Geometric view factors front->roof done (13/45)')

% From the right-hand-side cubicle wall (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(1)>pjjj(1) % visibility check            
            % Only if the roof patch is on the right of the cubicle wall,
            % there is visibility
            if piii(1)>(-1/2+2*d)
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
                            cosjjj = abs(tmp2(1));
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
    end
end
disp('Geometric view factors cubicleright->roof done (14/45)')

% From the left-hand-side cubicle wall (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the left-hand-side cubicle wall
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(1)<pjjj(1) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(n^2+iii,6*n^2+halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-left pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
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
                F(n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleleft->roof done (15/45)')

% From the cubicle top surface (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the cubicle top
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
                cosiii = abs(tmp2(2));
                cosjjj = abs(tmp2(3));
                % Evaluate integrand
                intgndmat(kkk,lll) = cosiii*cosjjj/(pi*r^2);
            end
        end
        % Calculate element of F
        viewfactor = qw*sum(sum(intgndmat))/d^2;
        F(n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;        
    end
end
disp('Geometric view factors cubicletop->roof done (16/45)')

% From the cubicle front surface (jjj) to the roof (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the roof
        piii = [Xmat(iii,2);Ymat(iii,2);Zmat(iii,2)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(2)<pjjj(2) % Visibility check
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
            F(n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubiclefront->roof done (17/45)')






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
        if abs(piii(1)-(-1/2+d))<d & piii(2)>0 % Floor pixels behind cubicle
            visible = 0;
        end
        % If the floor patch is to the left of the cubicle wall,
        % there is a chance of invisibility
        if piii(1)<(-1/2)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors right wall->floor done (18/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(2)>0 % Floor pixels behind cubicle
            visible = 0;
        end
        % Only if the floor patch is to the right of the cubicle wall,
        % there is a chance of invisibility
        if piii(1)>(-1/2+2*d)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors left->floor done (19/45)')

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
        if abs(piii(1)-(-1/2+d))<d & piii(2)>0 % Floor pixels behind cubicle
            visible = 0;
        end
        % Only if the two centerpoints are on opposite sides of cubicle wall,
        % there is a chance of invisibility
        if (piii(1)-(-1/2+d))*(pjjj(1)-(-1/2+d))<0
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors front->floor done (20/45)')

% From the right-hand-side cubicle wall (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(1)>pjjj(1) % visibility check            
            % Only if the floor patch is on the right of the cubicle wall,
            % there is visibility
            if piii(1)>(-1/2+2*d)
                % Distance between the points
                difvec0 = piii-pjjj;
                r0      = norm(difvec0);                
                % Check if the two pixels share an edge
                if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                    % Calculate element of F analytically
                    F(2*n^2+iii,6*n^2+jjj) = shared_edge_F;
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
                            % Quadrature point in the cubicle-right pixel
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
                    F(2*n^2+iii,6*n^2+jjj) = viewfactor;
                end
            end
        end
    end
end
disp('Geometric view factors cubicleright->floor done (21/45)')

% From the left-hand-side cubicle wall (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the left-hand-side cubicle wall
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(1)<pjjj(1) % visibility check
            % Distance between the points
            difvec0 = piii-pjjj;
            r0      = norm(difvec0);
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(2*n^2+iii,6*n^2+halfn^2+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-left pixel
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
                F(2*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleleft->floor done (22/45)')

% From the cubicle top surface (jjj) to the floor (iii)
% There is no visibility at all
disp('Geometric view factors cubicletop->floor (23/45)')

% From the cubicle front surface (jjj) to the floor (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the floor
        piii = [Xmat(iii,3);Ymat(iii,3);Zmat(iii,3)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        % Distance between the points
        difvec0 = piii-pjjj;
        r0      = norm(difvec0);
        if piii(2)<pjjj(2) % visibility check
            % Check if the two pixels share an edge
            if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                % Calculate element of F analytically
                F(2*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = shared_edge_F;
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
                        % Quadrature point in the cubicle-front pixel
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
                F(2*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubiclefront->floor done (24/45)')






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
        % Visibility check
        tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
        if ~(tmp3(2)>0 & tmp3(3)<0 )
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
end
disp('Geometric view factors left->right-hand-side wall done (25/45)')

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
disp('Geometric view factors front->right-hand-side wall done (26/45)')

% From the right-hand-side cubicle wall (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(1)>pjjj(1) % visibility check            
            % Only if the right-hand-side wall patch is on the right of the cubicle wall,
            % there is visibility
            if piii(1)>(-1/2+2*d)
                % Distance between the points
                difvec0 = piii-pjjj;
                r0      = norm(difvec0);                
                % Check if the two pixels share an edge
                if r0<(sqrt(2)*d/2 + 1e-8) % Edge shared
                    % Calculate element of F analytically
                    F(3*n^2+iii,6*n^2+jjj) = shared_edge_F;
                else % Edge not shared: integrate for F using quadrature
                    % Initialize matrix of integrand values at quadrature points
                    intgndmat = zeros(qn^2,qn^2);
                    % Double loop over four-dimensional quadrature
                    for kkk = 1:qn^2
                        for lll = 1:qn^2
                            % Quadrature point in the right-hand-side wall pixel
                            qpiii = piii;
                            qpiii(2) = qpiii(2)+q1(kkk);
                            qpiii(3) = qpiii(3)+q2(kkk);
                            % Quadrature point in the cubicle-right pixel
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
                    F(3*n^2+iii,6*n^2+jjj) = viewfactor;
                end
            end
        end
    end
end
disp('Geometric view factors cubicleright->right-hand-side wall done (27/45)')

% From the left-hand-side cubicle wall (jjj) to the right-hand-side wall (iii)
% There is no visibility at all
disp('Geometric view factors cubicleleft->right-hand-side wall done (28/45)')

% From the cubicle top surface (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the cubicle top
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
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
                    % Quadrature point in the cubicle-top pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
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
            F(3*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubicletop->right-hand-side wall (29/45)')

% From the cubicle front surface (jjj) to the right-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the right wall
        piii = [Xmat(iii,4);Ymat(iii,4);Zmat(iii,4)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        if piii(2)<pjjj(2) % visibility check
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
                    % Quadrature point in the cubicle-front pixel
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
            F(3*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubiclefront->right-hand-side wall done (30/45)')
        






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
        % Only if the front patch is on the right of the cubicle wall,
        % there is a chance of invisibility
        if pjjj(1)>(-1/2+2*d)
            tmp3 = pjjj+abs((-1/2+d-pjjj(1))/tmp2(1))*tmp2;
            if tmp3(2)>0 & tmp3(3)<0 % Does the connecting line cross the cubicle wall?
                visible = 0;
            end
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
disp('Geometric view factors front->left-hand-side wall done (31/45)')

% From the right-hand-side cubicle wall (jjj) to the left-hand-side wall (iii)
% No visibility
disp('Geometric view factors cubicleright->left-hand-side wall done (32/45)')

% From the left-hand-side cubicle wall (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in the left-hand-side cubicle wall
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(1)<pjjj(1) % visibility check
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
                        % Quadrature point in the cubicle-left pixel
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
                F(4*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleleft->left-hand-side wall done (33/45)')

% From the cubicle top surface (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in the cubicle top
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
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
                    % Quadrature point in the cubicle-top pixel
                    qpjjj = pjjj;
                    qpjjj(1) = qpjjj(1)+q1(lll);
                    qpjjj(2) = qpjjj(2)+q2(lll);
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
            F(4*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubicletop->left-hand-side wall (34/45)')

% From the cubicle front surface (jjj) to the left-hand-side wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the left wall
        piii = [Xmat(iii,5);Ymat(iii,5);Zmat(iii,5)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
        if piii(2)<pjjj(2) % visibility check
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
                    % Quadrature point in the cubicle-front pixel
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
            F(4*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubiclefront->left-hand-side wall done (35/45)')
        





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From the right-hand-side cubicle wall (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in the right-hand-side cubicle wall
        pjjj = [Xmat2(jjj,1);Ymat2(jjj,1);Zmat2(jjj,1)];
        if piii(1)>pjjj(1) % visibility check
            % Only if the back wall patch is on the right of the cubicle wall,
            % there is visibility
            if piii(1)>(-1/2+2*d)
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
                        % Quadrature point in the cubicle-right pixel
                        qpjjj = pjjj;
                        qpjjj(2) = qpjjj(2)+q1(lll);
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
                F(5*n^2+iii,6*n^2+jjj) = viewfactor;
            end
        end
    end
end
disp('Geometric view factors cubicleright->front wall done (36/45)')

% From the left-hand-side cubicle wall (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:halfn^2
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in the left-hand-side cubicle wall
        pjjj = [Xmat2(jjj,2);Ymat2(jjj,2);Zmat2(jjj,2)];
        if piii(1)<pjjj(1) % visibility check
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
                    % Quadrature point in the cubicle-left pixel
                    qpjjj = pjjj;
                    qpjjj(2) = qpjjj(2)+q1(lll);
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
            F(5*n^2+iii,6*n^2+halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubicleleft->front wall done (37/45)')

% From the cubicle top surface (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in the cubicle top
        pjjj = [Xmat3(jjj,1);Ymat3(jjj,1);Zmat3(jjj,1)];
        if piii(3)>pjjj(3) % visibility check
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
                    % Quadrature point in the cubicle-top pixel
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
            F(5*n^2+iii,6*n^2+2*halfn^2+jjj) = viewfactor;
        end
    end
end
disp('Geometric view factors cubicletop->front wall (38/45)')

% From the cubicle front surface (jjj) to the front wall (iii)
for iii = 1:n^2
    for jjj = 1:2*halfn
        % Centerpoint of the current pixel in the front wall
        piii = [Xmat(iii,6);Ymat(iii,6);Zmat(iii,6)];
        % Centerpoint of the current pixel in the cubicle front
        pjjj = [Xmat3(jjj,2);Ymat3(jjj,2);Zmat3(jjj,2)];
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
                % Quadrature point in the cubicle-front pixel
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
        F(5*n^2+iii,6*n^2+2*halfn^2+2*halfn+jjj) = viewfactor;
    end
end
disp('Geometric view factors cubiclefront->front wall done (39/45)')

% The sides of the cubicle wall have no visibility to each other, so 
% (40-45/45) are all zeros.
disp('The rest of geometric view factors (40-45/45)')


% Use symmetry to finish the construction of F. F is symmetric
% since all the pixels in our model have equal size.
F = F+F.';

% Check the matrix F. Most of the row sums should be one; this fact can be 
% used in debugging. 
figure(10)
clf
plot(sum(F))
% figure(11)
% clf
% spy(F)


% Save matrix to disc
save data/F_wall F n qn d Xmat Ymat Zmat Xmat2 Ymat2 Zmat2 Xmat3 Ymat3 Zmat3 halfn n_wall n_cubicle
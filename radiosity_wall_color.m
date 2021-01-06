% Further testing for building a radiosity lighting method for virtual
% spaces.
%
% Samuli Siltanen November 2020

%% Preliminaries

% Load precomputed stuff
disp('Loading data')
load data/F_wall F n qn d Xmat Ymat Zmat Xmat2 Ymat2 Zmat2 Xmat3 Ymat3 Zmat3 halfn n_wall n_cubicle
disp('Data loaded')

% Adjust the dark shades. Colors darker than the threshold will become
% black, so increasing the threshold will darken the image. 
threshold = 0.03;

% Sigmoid correction for optimal gray levels. Increasing betapar1 will
% darken the image, especially the shadows. Increasing betapar2 will
% lighten the image, especially highlights. 
betapar1 = 1.4;
betapar2 = 6;

% Camera settings. The camera is located at vector "campos", it is pointed
% towards "camtar", and the view angle is "camang". A larger "camang" value
% will give a more "wide-angle lens", so more of the scene is seen in the
% image. 
campos = [.2 -2.3 -.30];
camtar = [0 0 0];
camang = 70;

%% Construct the color vector (B-vector) using the radiosity lighting model.

% Construct the right hand side Evec of the radiosity equation. Evec
% describes the contribution of emitted light in the scene. For example,
% each pixel belonging to a lamp in the virtual space causes a positive
% element in Evec.
Evec   = zeros(n_wall+n_cubicle,1);
indvec = repmat(logical(0),size(Evec));
indvec(n^2+[1:n^2]) = sqrt((Xmat(:,2)-.3).^2+Ymat(:,2).^2)<.3; % Ceiling lamp
indvec(4*n^2+[1:n^2]) = ...
    ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)-0)<.1))|...
    ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)-1/2)<.1))|...
    ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)+1/2)<.1)); % Rectangular lamps in the left wall
% indvec(4*n^2+[1:n^2]) = ...
%     sqrt((Zmat(:,5)-0).^2+(Ymat(:,5)-0).^2)<.1 |...
%     sqrt((Zmat(:,5)-0).^2+(Ymat(:,5)-1/2).^2)<.1 |...
%     sqrt((Zmat(:,5)-0).^2+(Ymat(:,5)+1/2).^2)<.1; % Round lamps in the left wall
% Evec(n^2+round(n^2/2)-2) = 1;
% Evec(3*n^2+round(n^2/2)-2) = 1;
Evec(indvec) = 1;
disp('Right-hand-side constructed')

% The parameter rho adjusts the surface material (how much incoming light
% is reflected away from a patch, 0<rho<=1)
rho = .9*ones(n_wall+n_cubicle,1);
rho(n^2+[1:n^2]) = 1; % Bright ceiling
rho(2*n^2+[1:n^2]) = .7; % Dark floor

% Solve for color vector.
disp('Solving radiosity equation...')
tic
colorvec_orig = gmres(eye(n_wall+n_cubicle)-repmat(rho,1,n_wall+n_cubicle).*F,Evec);
disp(['Radiosity equation solved in ',num2str(toc),' seconds'])


%% Produce a still image of the scene, for checking

% Adjust the dark shades and normalize the values of the color vector
% between 0 and 1.
colorvec = colorvec_orig-threshold;
colorvec = max(0,colorvec);
colorvec = colorvec/max(colorvec);

% Sigmoid correction for optimal gray levels.
colorvec  = betacdf(colorvec,betapar1,betapar2);
% figure(100)
% clf
% t = linspace(0,1,200);
% plot(t,betacdf(t,betapar1,betapar2));
% title('Grayscale adjustment curve')
% Construct grayscale color matrix by repeating the same color vector for
% red, green and blue channels.


% Construct color matrix
colorvecR = colorvec;
colorvecR(1:n^2) = colorvecR(1:n^2).^.8; % Back wall
colorvecR(3*n^2+[1:n^2]) = colorvecR(3*n^2+[1:n^2]).^1.2; % Right wall
colorvecR(4*n^2+[1:n^2]) = colorvecR(4*n^2+[1:n^2]).^1.2; % Left wall
colorvecG = colorvec;
colorvecG(1:n^2) = colorvecG(1:n^2).^1.4; % Back wall
colorvecB = colorvec;
colorvecB(1:n^2) = colorvecB(1:n^2).^1.4; % Back wall
colorvecB(3*n^2+[1:n^2]) = colorvecB(3*n^2+[1:n^2]).^.7; % Right wall
colorvecB(4*n^2+[1:n^2]) = colorvecB(4*n^2+[1:n^2]).^.7; % Left wall
colormat = [colorvecR(:),colorvecG(:),colorvecB(:)];



% Create plot window
figure(1)
clf

% Draw all the walls consisting of n x n little squares (pixels).
% Pick the gray value of each square from the illumination vector
% calculated by the radiosity method above

% The back wall
colorind = 1;
for iii = 1:(n^2)
    p1 = patch(...
        [Xmat(iii,1)+d/2,Xmat(iii,1)+d/2,Xmat(iii,1)-d/2,Xmat(iii,1)-d/2],...
        [Ymat(iii,1),Ymat(iii,1),Ymat(iii,1),Ymat(iii,1)],...
        [Zmat(iii,1)-d/2,Zmat(iii,1)+d/2,Zmat(iii,1)+d/2,Zmat(iii,1)-d/2],...
        colormat(colorind,:));
    hold on
    % If the following line is missing, each patch will have black
    % boundaries. Here we set the boundary color to match the patch color,
    % so no boundaries are visible.
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Roof
for iii = 1:(n^2)
    p1 = patch(...
        [Xmat(iii,2)+d/2,Xmat(iii,2)+d/2,Xmat(iii,2)-d/2,Xmat(iii,2)-d/2],...
        [Ymat(iii,2)-d/2,Ymat(iii,2)+d/2,Ymat(iii,2)+d/2,Ymat(iii,2)-d/2],...
        [Zmat(iii,2),Zmat(iii,2),Zmat(iii,2),Zmat(iii,2)],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Floor
for iii = 1:(n^2)
    p1 = patch(...
        [Xmat(iii,3)+d/2,Xmat(iii,3)+d/2,Xmat(iii,3)-d/2,Xmat(iii,3)-d/2],...
        [Ymat(iii,3)-d/2,Ymat(iii,3)+d/2,Ymat(iii,3)+d/2,Ymat(iii,3)-d/2],...
        [Zmat(iii,3),Zmat(iii,3),Zmat(iii,3),Zmat(iii,3)],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Right-hand-side wall
for iii = 1:(n^2)
    p1 = patch(...
        [Xmat(iii,4),Xmat(iii,4),Xmat(iii,4),Xmat(iii,4)],...
        [Ymat(iii,4)+d/2,Ymat(iii,4)+d/2,Ymat(iii,4)-d/2,Ymat(iii,4)-d/2],...
        [Zmat(iii,4)-d/2,Zmat(iii,4)+d/2,Zmat(iii,4)+d/2,Zmat(iii,4)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Left-hand-side wall
for iii = 1:(n^2)
    p1 = patch(...
        [Xmat(iii,5),Xmat(iii,5),Xmat(iii,5),Xmat(iii,5)],...
        [Ymat(iii,5)+d/2,Ymat(iii,5)+d/2,Ymat(iii,5)-d/2,Ymat(iii,5)-d/2],...
        [Zmat(iii,5)-d/2,Zmat(iii,5)+d/2,Zmat(iii,5)+d/2,Zmat(iii,5)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Front wall (not drawn)
for iii = 1:(n^2)
    %     p1 = patch(...
    %         [Xmat(iii,6),Xmat(iii,6),Xmat(iii,6),Xmat(iii,6)],...
    %         [Ymat(iii,6)+d/2,Ymat(iii,6)+d/2,Ymat(iii,6)-d/2,Ymat(iii,6)-d/2],...
    %         [Zmat(iii,6)-d/2,Zmat(iii,6)+d/2,Zmat(iii,6)+d/2,Zmat(iii,6)-d/2],...
    %         colormat(colorind,:));
    %     set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Right-hand-side CUBICLE wall
for iii = 1:(halfn^2)
    p1 = patch(...
        [Xmat2(iii,1),Xmat2(iii,1),Xmat2(iii,1),Xmat2(iii,1)],...
        [Ymat2(iii,1)+d/2,Ymat2(iii,1)+d/2,Ymat2(iii,1)-d/2,Ymat2(iii,1)-d/2],...
        [Zmat2(iii,1)-d/2,Zmat2(iii,1)+d/2,Zmat2(iii,1)+d/2,Zmat2(iii,1)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Left-hand-side CUBICLE wall
for iii = 1:(halfn^2)
    p1 = patch(...
        [Xmat2(iii,2),Xmat2(iii,2),Xmat2(iii,2),Xmat2(iii,2)],...
        [Ymat2(iii,2)+d/2,Ymat2(iii,2)+d/2,Ymat2(iii,2)-d/2,Ymat2(iii,2)-d/2],...
        [Zmat2(iii,2)-d/2,Zmat2(iii,2)+d/2,Zmat2(iii,2)+d/2,Zmat2(iii,2)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Top surface of CUBICLE wall
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,1)-d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)-d/2],...
        [Ymat3(iii,1)+d/2,Ymat3(iii,1)+d/2,Ymat3(iii,1)-d/2,Ymat3(iii,1)-d/2],...
        [Zmat3(iii,1),Zmat3(iii,1),Zmat3(iii,1),Zmat3(iii,1)],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Front surface of CUBICLE wall
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,2)-d/2,Xmat3(iii,2)+d/2,Xmat3(iii,2)+d/2,Xmat3(iii,2)-d/2],...
        [Ymat3(iii,2),Ymat3(iii,2),Ymat3(iii,2),Ymat3(iii,2)],...
        [Zmat3(iii,2)+d/2,Zmat3(iii,2)+d/2,Zmat3(iii,2)-d/2,Zmat3(iii,2)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Camera settings
camproj('perspective')
set(gca,'CameraPosition',campos,'CameraTarget',camtar,'CameraViewAngle',camang)

% Axis settings
axis equal
axis off
drawnow









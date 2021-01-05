% Radiosity lighting method for a virtual room, in color.
% The routine "radiosity_emptyroom_Fcomp.m" needs to be computed before
% this one. 
%
% Samuli Siltanen January 2021

%% Preliminaries

% Load precomputed stuff
disp('Loading data')
load data/F_emptyroom F n qn d Xmat Ymat Zmat 
disp('Data loaded')



%% Construct the color vector (B-vector) using the radiosity lighting model.

% Construct the right hand side Evec of the radiosity equation. Evec
% describes the contribution of emitted light in the scene. For example,
% each pixel belonging to a lamp in the virtual space causes a positive
% element in Evec.
Evec   = zeros(6*n^2,1);
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
rho = .9*ones(6*n^2,1);
rho(n^2+[1:n^2]) = 1; % Bright ceiling
rho(2*n^2+[1:n^2]) = .7; % Dark floor

% Solve for color vector.
disp('Solving radiosity equation...')
tic
colorvec_orig = gmres(eye(6*n^2)-repmat(rho,1,6*n^2).*F,Evec);
disp(['Radiosity equation solved in ',num2str(toc),' seconds'])


%% Produce a still image of the scene

% Adjust the dark shades and normalize the values of the color vector 
% between 0 and 1.
threshold = 0.05;
colorvec = colorvec_orig-threshold;
colorvec = max(0,colorvec);
colorvec = colorvec/max(colorvec);

% Sigmoid correction for optimal gray levels.
betapar1 = 1;
betapar2 = 10;
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

% Camera settings
camproj('perspective')
set(gca,'CameraPosition',[-.2 -3 -.30],'CameraTarget',[0 0 .1],'CameraViewAngle',50)

% Axis settings
axis equal
axis off
drawnow







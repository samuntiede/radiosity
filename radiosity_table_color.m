% Further testing for building a radiosity lighting method for virtual
% spaces. Here the room has a table. 
%
% Samuli Siltanen January 2021

%% Preliminaries

% Load precomputed stuff
disp('Loading data')
load data/F_table F n qn d Xmat Ymat Zmat Xmat2 Ymat2 Zmat2 Xmat3 Ymat3 Zmat3 halfn n_wall n_table t_top t_front t_left
disp('Data loaded')

% Adjust the dark shades. Colors darker than the threshold will become
% black, so increasing the threshold will darken the image.
threshold = 0.02;

% Sigmoid correction for optimal gray levels. Increasing betapar1 will
% darken the image, especially the shadows. Increasing betapar2 will
% lighten the image, especially highlights. 
betapar1 = 2.2;
betapar2 = 10;

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
Evec   = zeros(n_wall+n_table,1);
indvec = repmat(logical(0),size(Evec));
indvec(n^2+[1:n^2]) = sqrt((Xmat(:,2)-.3).^2+Ymat(:,2).^2)<.3; % Ceiling lamp
% indvec(4*n^2+[1:n^2]) = ...
%     ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)-0)<.1))|...
%     ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)-1/2)<.1))|...
%     ((abs(Zmat(:,5)-0)<.3)&(abs(Ymat(:,5)+1/2)<.1)); % Rectangular lamps in the left wall
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
rho = .9*ones(n_wall+n_table,1);
rho(n^2+[1:n^2]) = 1; % Bright ceiling
rho(2*n^2+[1:n^2]) = .7; % Dark floor

% Solve for color vector.
disp('Solving radiosity equation...')
tic
colorvec_orig = gmres(eye(n_wall+n_table)-repmat(rho,1,n_wall+n_table).*F,Evec);
disp(['Radiosity equation solved in ',num2str(toc),' seconds'])

% Adjust the dark shades and normalize the values of the color vector
% between 0 and 1.
colorvec = colorvec_orig-threshold;
colorvec = max(0,colorvec);
colorvec = colorvec/max(colorvec);

% Sigmoid correction for optimal gray levels.
colorvec  = Scaled_BetaCDF(colorvec,betapar1,betapar2);
% figure(100)
% clf
% t = linspace(0,1,200);
% plot(t,Scaled_BetaCDF(t,betapar1,betapar2));
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




% Sigmoid correction for optimal gray levels.
betapar1 = 2;
betapar2 = 13;
colorvec  = Scaled_BetaCDF(colorvec,betapar1,betapar2);
% figure(100)
% clf
% t = linspace(0,1,200);
% plot(t,Scaled_BetaCDF(t,betapar1,betapar2));
% title('Grayscale adjustment curve')
% Construct grayscale color matrix by repeating the same color vector for
% red, green and blue channels.

% Random colors, for testing
% colorvec = .6*rand(n_wall+n_table,1);
% colorvec(6*n^2+[1:n_table]) = .5+.5*rand(n_table,1);

% Visibility colors, for testing
% colorvec = F(:,2.5*n^2+5);
% colorvec = colorvec/max(colorvec(:));
% colorvec = colorvec.^.2;

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

% Draw a background patch of "green-screen green" color
screengreen = [0 226 2]/255;
GSmax = 40;
GSSE = [GSmax;2;-GSmax];  % Front wall, southeast corner
GSSW = [-GSmax;2;-GSmax]; % Front wall, southwest corner
GSNW = [-GSmax;2;GSmax];  % Front wall, northwest corner
GSNE = [GSmax;2;GSmax];   % Front wall, northeast corner
%patch([GSSE(1);GSSW(1);GSNW(1);GSNE(1)],[GSSE(2);GSSW(2);GSNW(2);GSNE(2)],[GSSE(3);GSSW(3);GSNW(3);GSNE(3)], screengreen)
hold on

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


% Table top
for iii = 1:(halfn^2)
    p1 = patch(...
        [Xmat2(iii,1)-d/2,Xmat2(iii,1)+d/2,Xmat2(iii,1)+d/2,Xmat2(iii,1)-d/2],...
        [Ymat2(iii,1)+d/2,Ymat2(iii,1)+d/2,Ymat2(iii,1)-d/2,Ymat2(iii,1)-d/2],...
        [Zmat2(iii,1),Zmat2(iii,1),Zmat2(iii,1),Zmat2(iii,1)],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Table bottom
for iii = 1:(halfn^2)
    p1 = patch(...
        [Xmat2(iii,2)-d/2,Xmat2(iii,2)+d/2,Xmat2(iii,2)+d/2,Xmat2(iii,2)-d/2],...
        [Ymat2(iii,2)+d/2,Ymat2(iii,2)+d/2,Ymat2(iii,2)-d/2,Ymat2(iii,2)-d/2],...
        [Zmat2(iii,2),Zmat2(iii,2),Zmat2(iii,2),Zmat2(iii,2)],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Back side of table
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,1)-d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)-d/2],...
        [Ymat3(iii,1),Ymat3(iii,1),Ymat3(iii,1),Ymat3(iii,1)],...
        [Zmat3(iii,1)+d/2,Zmat3(iii,1)+d/2,Zmat3(iii,1)-d/2,Zmat3(iii,1)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Right side of table
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,2),Xmat3(iii,2),Xmat3(iii,2),Xmat3(iii,2)],...
        [Ymat3(iii,2)-d/2,Ymat3(iii,2)+d/2,Ymat3(iii,2)+d/2,Ymat3(iii,2)-d/2],...
        [Zmat3(iii,2)+d/2,Zmat3(iii,2)+d/2,Zmat3(iii,2)-d/2,Zmat3(iii,2)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Left side of table
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,3),Xmat3(iii,3),Xmat3(iii,3),Xmat3(iii,3)],...
        [Ymat3(iii,3)-d/2,Ymat3(iii,3)+d/2,Ymat3(iii,3)+d/2,Ymat3(iii,3)-d/2],...
        [Zmat3(iii,3)+d/2,Zmat3(iii,3)+d/2,Zmat3(iii,3)-d/2,Zmat3(iii,3)-d/2],...
        colormat(colorind,:));
    set(p1,'EdgeColor',colormat(colorind,:))
    colorind = colorind+1;
end

% Front side of table
for iii = 1:(2*halfn)
    p1 = patch(...
        [Xmat3(iii,4)-d/2,Xmat3(iii,4)+d/2,Xmat3(iii,4)+d/2,Xmat3(iii,4)-d/2],...
        [Ymat3(iii,4),Ymat3(iii,4),Ymat3(iii,4),Ymat3(iii,4)],...
        [Zmat3(iii,4)+d/2,Zmat3(iii,4)+d/2,Zmat3(iii,4)-d/2,Zmat3(iii,4)-d/2],...
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



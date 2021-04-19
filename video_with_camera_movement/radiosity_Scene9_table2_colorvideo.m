% Further testing for building a radiosity lighting method for virtual
% spaces.
%
% Samuli Siltanen November 2020

%% Preliminaries
Nframes = 300; % Must be even
t = linspace(0,1,Nframes/2);
% t = t.^(1.2); % Ease
% t = t(:);
% tt = pi*[t-1;1-flipud(t)];
% camposxvec = .2*cos(tt);
% camposzvec = .15*sin(tt);
t = t(:).^(2); % Ease
tt = [t;2-flipud(t)]/2;
camposxvec = -.1+tt*(.1-(-.1));
camposyvec = -3+tt*.3;
camposzvec = -.2+tt*(-.1);
camtarxvec = tt*(.1);
camtaryvec = tt*.1;
camtarzvec = .15+tt*(.02);
camviewanglevec = 50+tt*10;

% Background patch of "green-screen green" color
screengreen = [0 226 2]/255;
GSmax = 40;
GSSE = [GSmax;2;-GSmax];  % Front wall, southeast corner
GSSW = [-GSmax;2;-GSmax]; % Front wall, southwest corner
GSNW = [-GSmax;2;GSmax];  % Front wall, northwest corner
GSNE = [GSmax;2;GSmax];   % Front wall, northeast corner

% Load precomputed stuff
disp('Loading data')
load data/F_table2 F n qn d Xmat Ymat Zmat Xmat2 Ymat2 Zmat2 Xmat3 Ymat3 Zmat3 halfn n_wall n_table t_top t_front t_left
disp('Data loaded')

% Adjust the dark shades
threshold = 0.02;

% Sigmoid correction for optimal gray levels.
betapar1 = 2;
betapar2 = 15;

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



%% Create video



% Open video file
videofilename = 'radiosity_Scene9_table2_video2_color';
videotype = 'MPEG-4';
v1 = VideoWriter(videofilename,videotype);
v1.Quality = 95;
open(v1);


% Create plot window
figure(1)
clf

% Draw green screen
patch([GSSE(1);GSSW(1);GSNW(1);GSNE(1)],[GSSE(2);GSSW(2);GSNW(2);GSNE(2)],[GSSE(3);GSSW(3);GSNW(3);GSNE(3)], screengreen)
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

% 
% % Table top
% for iii = 1:(halfn^2)
%     p1 = patch(...
%         [Xmat2(iii,1)-d/2,Xmat2(iii,1)+d/2,Xmat2(iii,1)+d/2,Xmat2(iii,1)-d/2],...
%         [Ymat2(iii,1)+d/2,Ymat2(iii,1)+d/2,Ymat2(iii,1)-d/2,Ymat2(iii,1)-d/2],...
%         [Zmat2(iii,1),Zmat2(iii,1),Zmat2(iii,1),Zmat2(iii,1)],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end
% 
% % Table bottom
% for iii = 1:(halfn^2)
%     p1 = patch(...
%         [Xmat2(iii,2)-d/2,Xmat2(iii,2)+d/2,Xmat2(iii,2)+d/2,Xmat2(iii,2)-d/2],...
%         [Ymat2(iii,2)+d/2,Ymat2(iii,2)+d/2,Ymat2(iii,2)-d/2,Ymat2(iii,2)-d/2],...
%         [Zmat2(iii,2),Zmat2(iii,2),Zmat2(iii,2),Zmat2(iii,2)],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end
% 
% % Back side of table
% for iii = 1:(2*halfn)
%     p1 = patch(...
%         [Xmat3(iii,1)-d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)+d/2,Xmat3(iii,1)-d/2],...
%         [Ymat3(iii,1),Ymat3(iii,1),Ymat3(iii,1),Ymat3(iii,1)],...
%         [Zmat3(iii,1)+d/2,Zmat3(iii,1)+d/2,Zmat3(iii,1)-d/2,Zmat3(iii,1)-d/2],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end
% 
% % Right side of table
% for iii = 1:(2*halfn)
%     p1 = patch(...
%         [Xmat3(iii,2),Xmat3(iii,2),Xmat3(iii,2),Xmat3(iii,2)],...
%         [Ymat3(iii,2)-d/2,Ymat3(iii,2)+d/2,Ymat3(iii,2)+d/2,Ymat3(iii,2)-d/2],...
%         [Zmat3(iii,2)+d/2,Zmat3(iii,2)+d/2,Zmat3(iii,2)-d/2,Zmat3(iii,2)-d/2],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end
% 
% % Left side of table
% for iii = 1:(2*halfn)
%     p1 = patch(...
%         [Xmat3(iii,3),Xmat3(iii,3),Xmat3(iii,3),Xmat3(iii,3)],...
%         [Ymat3(iii,3)-d/2,Ymat3(iii,3)+d/2,Ymat3(iii,3)+d/2,Ymat3(iii,3)-d/2],...
%         [Zmat3(iii,3)+d/2,Zmat3(iii,3)+d/2,Zmat3(iii,3)-d/2,Zmat3(iii,3)-d/2],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end
% 
% % Front side of table
% for iii = 1:(2*halfn)
%     p1 = patch(...
%         [Xmat3(iii,4)-d/2,Xmat3(iii,4)+d/2,Xmat3(iii,4)+d/2,Xmat3(iii,4)-d/2],...
%         [Ymat3(iii,4),Ymat3(iii,4),Ymat3(iii,4),Ymat3(iii,4)],...
%         [Zmat3(iii,4)+d/2,Zmat3(iii,4)+d/2,Zmat3(iii,4)-d/2,Zmat3(iii,4)-d/2],...
%         colormat(colorind,:));
%     set(p1,'EdgeColor',colormat(colorind,:))
%     colorind = colorind+1;
% end

% Loop over frames
disp('Loop over frames')

for fff = 1:Nframes
    
    
    % Camera settings
    camproj('perspective')
    set(gca,'CameraPosition',[camposxvec(fff) camposyvec(fff) camposzvec(fff)],...
        'CameraTarget',[camtarxvec(fff) camtaryvec(fff) camtarzvec(fff)],...
        'CameraViewAngle',camviewanglevec(fff))
    
    % Axis settings
    axis equal
    axis off
    drawnow
    
    % Initial image
    im1 = print('-r300','-RGBImage');
    [row1,col1] = size(im1(:,:,1));
    
    % Adjust image height to 1080
    im2 = uint8(ones(1080,1920,3));
    im2(:,:,1) = uint8(255*screengreen(1));
    im2(:,:,2) = uint8(255*screengreen(2));
    im2(:,:,3) = uint8(255*screengreen(3));
    im3 = imresize(im1, [1080 NaN]);
    [row3, col3, tmp3] = size(im3);
    im2(:,round((1920-col3)/2)+[1:col3],1) = im3(:,:,1);
    im2(:,round((1920-col3)/2)+[1:col3],2) = im3(:,:,2);
    im2(:,round((1920-col3)/2)+[1:col3],3) = im3(:,:,3);
    
    
    % Set lamp pixels to white. For some reason the lamp appears as perfect
    % black
    im2R = im2(:,:,1);
    im2G = im2(:,:,2);
    im2B = im2(:,:,3);
    lamp_ind = im2R==0 & im2G==0 & im2B==0 ;
    im2R(lamp_ind) = 255;
    im2G(lamp_ind) = 255;
    im2B(lamp_ind) = 255;
    im2(:,:,1) = im2R;
    im2(:,:,2) = im2G;
    im2(:,:,3) = im2B;
    
    %     figure(3)
    %     clf
    %     imshow(im2)
    %     drawnow
    
    % Add frame to video
    writeVideo(v1,im2);
    
    % Monitor the run
    disp([fff Nframes])
end

close(v1);

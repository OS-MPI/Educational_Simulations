function [ImageOut] = X_Space2D_Demo(ImageIn,PlottingOn,SaveVid,varargin)
%% Introduction

% This function is a demonstration of the working principles behind X-Space
% Magnetic Particle Imaging (MPI). It makes a number of assumptions to
% simplify the system to make it easier to read and therefore learn from.
% Namely, the magnetic fields are assumed to be ideal (drive is
% homogeneous, and the FFP field is an ideal linear gradient).
% This function was written as a teaching tool and you should feel free to
% use it as such.

% The inputs are:
%   ImageIn - A path to the object you want to simulate, or the object itself. It should be
%   a 2D matrix of numeric values
%   PlottingOn - Logical (1=on, 0 =off) that either enables or disables the
%   plotting of image as it is created.
%   SaveVid - if it is 1, the function will save a video of the plotting
%   with path of "ImageIn_2DXSpaceVid" if the function is given an image
%   path in, or default to "[Date]_2DXSpaceVid"
addpath([pwd,'\Called Functions'])

NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments in (past "SaveVid") must be even')
elseif NumVarArgsIn==0
        
        AddNoiseSelected = 0;
        FilterFundamentalSelect = 0;
        PlotRateSelect = 0;
end
varargin = reshape(varargin,2,NumVarArgsIn/2);

for i = 1:NumVarArgsIn/2
   if strcmp(varargin(1,i),'AddNoise')
       AddNoise = cell2mat(varargin(2,i));
       AddNoiseSelected = 1; %Indicate the user selected to add noise or not
   elseif exist('AddNoiseSelected','var')
   else
        AddNoiseSelected = 0;
   end
   
   if strcmp(varargin(1,i),'FilterFundamental')
       FilterFundamental = cell2mat(varargin(2,i));
       FilterFundamentalSelect = 1; %Indicate the user selected to add noise or not
   elseif exist('FilterFundamentalSelect','var')
   else
        FilterFundamentalSelect = 0;
   end
   if strcmp(varargin(1,i),'PlotRate')
       PlotRate = cell2mat(varargin(2,i));
       PlotRateSelect = 1; %Indicate the user selected to add noise or not
   elseif exist('PlotRateSelect','var')
   else
        PlotRateSelect = 0;
   end   
end
if isstring(ImageIn) || ischar(ImageIn)
    ImagePath = ImageIn;
    ExtensionStartIndex = find(ImagePath=='.');
    ImagePath(ExtensionStartIndex:end) = []; %Cropping the file extension for later use
    ImageIn = imread(ImageIn);
    if ndims(ImageIn)==3
        ImageIn = sum(abs(ImageIn),3);
    end
    
else
    ImagePath = date();
    ImagePath(ImagePath=='-')='_';
end
% ImageIn = ImageIn>max(max(ImageIn))/2; % If you want to make the object binary uncomment this

fDriveX = 25.5e3;% Drive Frequency, Hz

Fs= 1e6; % Sampling rate, Hz
PointsPerDrivePrd = round(Fs/fDriveX);
x = -.03:.001:.03;%FOV meters

DrivePeriods = length(x)*.5;
% y = -.005:.0005:.005;%FOV meters
y = x;
[X,Y] = meshgrid(x,y);
ConcentrationMap = imresize(ImageIn,size(X));
ConcentrationMap = cast(ConcentrationMap,'double');
Grad = 1; %Gradient in Tesla per meter, this is somewhat arbitrary as the particle simuation isn't very accurate
BfflFunc = @(X) Grad.*X; %Set to be an ideal linear gradient, this can be modified to make the gradient non-ideal
B_FFP_X_TMP = BfflFunc(X); %B field of the FFP
B_FFP_Y_TMP = BfflFunc(Y); %B field of the FFP

t = 0:1/Fs:DrivePeriods*1/fDriveX; %Creation of the time vector
dt = t(2)-t(1);

tSig = t(1)+dt/2:dt:t(end)-dt/2; %The signal vector is offset in time by dt/2 due to the differentiation, so I make a new time vector to correspond to it


% BDriveX(1,1,:) = sawtooth(2*pi*fDriveX*t,.5)/2+sawtooth(2*pi*1/(t(end)+dt)*t); % A more complex drive field pattern
BDriveX(1,1,:) = sin(2*pi*fDriveX*t);
BDriveX = repmat(BDriveX,length(y),length(x));


BDriveY(1,1,:) = sawtooth(2*pi*1/(t(end)+dt)*t);
% BDriveY(1,1,:) = sawtooth(2*pi*3/(t(end)+dt)*t,.5); %A more complex Shift field pattern
BDriveY = repmat(BDriveY,length(y),length(x));



B_FFP_X = repmat(B_FFP_X_TMP,1,1,length(t));
BX = B_FFP_X+(1.1*Grad*x(end)*BDriveX); %This line defines the B field for all points in time in the X direction,
%it is the summation of the FFP field plus the drive field. the 1.1 makes sure the FFP covers the full field of view
B_FFP_Y = repmat(B_FFP_Y_TMP,1,1,length(t));
BY = B_FFP_Y+(1.1*Grad*y(end)*BDriveY);

clear B_FFP_Y_TMP B_FFP_X_TMP


BMag = sqrt(BX.^2+BY.^2); %Magnitude of B Field
FFP_Loc_Orig_X = X(BMag(:,:,1)==min(min(BMag(:,:,1)))); %Original locations of the FFP is the minimum of the magnitude
FFP_Loc_Orig_Y = Y(BMag(:,:,1)==min(min(BMag(:,:,1))));


 Mx_Demo = squeeze(sum(sum(Langevin(BMag,ConcentrationMap).*BX./(BMag),2),1));
    
    
    Vx_Demo = diff(BX(1,1,:)); %Taking the velocity at an arbitrary location
    Vx_Demo=Vx_Demo(:); %making a column

    Pos_X = FFP_Loc_Orig_X-cumtrapz(Vx_Demo/Grad);

    
    Vy_Demo = diff(BY(1,1,:));
    Vy_Demo = Vy_Demo(:);
    Pos_Y = -FFP_Loc_Orig_Y+cumtrapz(Vy_Demo/Grad);
    
    SigX_Demo = diff(Mx_Demo(:));
    CroppingFactor = 0.3; %The cropping factor controls how much data is deleted due to being at low velocities. 0.3= 30% of data is deleted. 
    HighVelocityIndex = abs(Vx_Demo)>(CroppingFactor*max(Vx_Demo(:))); %This is just a vector with ones corresponding to the times where the velocity is above 0.3* the max
    

    SigX_Demo(not(HighVelocityIndex)) = NaN;
    
    
    MaxSig = max(SigX_Demo(:)./Vx_Demo(:)); %Vel Cor
    %         MaxSig = max(abs(SigX_Demo(:))); %not Vel Cor
    
    MinSig = min(SigX_Demo(:)./Vx_Demo(:)); %Vel Cor
    %         MinSig = min(abs(SigX_Demo(:))); %Not Vel COr




CometOn = 1; %An option for plotting, must be 0 or 1

if PlottingOn==1
    if SaveVid==1
        VideoPath = [ImagePath,'_2DXSpaceVid'];
    end
    
    figure('Position',[100 100 700 700])
    
    Ax2 = axes('Position',[0.1 0.725 0.8 0.25]);
    % Ax4 = axes('Position',[0.05 0.62 0.9 0.125]);
    if CometOn==0
        Ax1 = axes('Position',[0.2 0.02 0.6 0.6]);
    elseif CometOn==1
        Ax1 = axes('Position',[0.05 0.15 0.4 0.4]);
        Ax5 = axes('Position',[0.55 0.15 0.4 0.4]);
    else
        error('Invalid CometOn Value')
    end
    
   
    
    
    %     EndPt = length(t);
    EndPt = length(t)-1;
    if PlotRateSelect==0
        PlotRate = 3; %This line controls how quickly the plotting occurs by skipping frames. If this is = 1, every frame is recorded, if it is equal to 10, every 10th frame is recorded
    end
    TimeVecDrivePrds = 10;
    TimeVecNumPoints = TimeVecDrivePrds*PointsPerDrivePrd;
    
        
    %%
    VidCount=1;
    for jj = 1:PlotRate:EndPt
        
        
        axes(Ax1)
        cla
        imagesc(sqrt(BX(:,:,jj).^2+BY(:,:,jj).^2))
        hold on
        
        quiver(BX(:,:,jj),BY(:,:,jj))
        title('Magnitude of B Field')
        caxis([0 max(max(BMag(:,:,1)))])
        
        set(gca,'XTick',[], 'YTick', [])
        colormap parula
        % legend('Nanoparticles')
        axes(Ax5)
        cla
        plot3(squeeze(Pos_X(1:jj)),Pos_Y(1:jj),squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj))) %Velocity Corrected
        
        scatter3(squeeze(Pos_X(1:jj)),Pos_Y(1:jj),squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj)),30,squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj)),'Filled') %Velocity Corrected
        
        % plot3(squeeze(Pos_X(1:jj)),Pos_Y(1:jj),abs(squeeze(SigX_Demo(1:jj)))) %Not Vel Correct
        
        hold on
        plot3(Pos_X(jj),Pos_Y(jj),squeeze(SigX_Demo(jj)./Vx_Demo(jj)),'r*') %Vel Cor
        % plot3(Pos_X(jj),Pos_Y(jj),squeeze(abs(SigX_Demo(jj))),'r*') %Not Vel Cor
        xlim([min(Pos_X) max(Pos_X)])
        ylim([min(Pos_Y) max(Pos_Y)])
        zlim([0 MaxSig])
        title({'Signal vs Position';'Vel. Corrected'})
        xlabel('X position')
        ylabel('Y Position')
        zlabel('VC Signal')
        
        axes(Ax2)
        cla
        plot(t(1:end-1)*1000,squeeze(SigX_Demo./Vx_Demo(:)),'b','LineWidth',1)
        hold on
        plot(t(jj)*1000,squeeze(SigX_Demo(jj)./Vx_Demo(jj)),'bo','LineWidth',2)
        
        
%         set(gca,'XTick',[])
        ylabel('Voltage')
        if jj<=ceil(TimeVecNumPoints/2)+1
        xlim([0 t(TimeVecNumPoints)*1000])
        elseif jj>ceil(TimeVecNumPoints/2) && jj<EndPt-floor(TimeVecNumPoints/2)
            xlim([t(jj-floor(TimeVecNumPoints/2)) t(jj+floor(TimeVecNumPoints/2))]*1000)
        else
            xlim([t(end-floor(TimeVecNumPoints)) t(end)]*1000)
        end
        ylim([MinSig MaxSig])
        xlabel('Time (milliseconds)','FontSize',12,'FontWeight','bold')
        xlabel('Time')
        
        
        
        
        pause(.001)
        if SaveVid==1
            frame_1(VidCount) = getframe(gcf);
            
            VidCount = VidCount+1;
        end
        
    end
    
    if SaveVid==1
        v = VideoWriter(VideoPath,'MPEG-4');
        open(v)
        writeVideo(v,frame_1)
        close(v)
    end
end
%%
if AddNoiseSelected == 0
    AddNoise = 0; %if this is 1, then Noise will be added at a magnitude determined by NoiseMag By default it is 0
end
NoiseMag = 0.5; %Magnitude of noise relative to the max signal 
NoisePct = NoiseMag*100;

if AddNoise==1
    SigX_Demo = SigX_Demo+rand(size(SigX_Demo))*NoiseMag*max(SigX_Demo(:));
end


figure,
jj = length(t)-1;

plot3(squeeze(Pos_X(1:jj)),Pos_Y(1:jj),squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj))) %Velocity COrrected
hold on
scatter3(squeeze(Pos_X(1:jj)),Pos_Y(1:jj),squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj)),30,squeeze(SigX_Demo(1:jj))./squeeze(Vx_Demo(1:jj)),'Filled') %Velocity Corrected
      
xlim([min(Pos_X) max(Pos_X)])
ylim([min(Pos_Y) max(Pos_Y)])
% zlim([0 MaxSig])
if AddNoise==1
    title({'Signal across space';['Vel. corrected w/',num2str(NoisePct),'% noise']})
else
    title({'Signal across space';['Vel. corrected w/ no added noise']})
end

xlabel('X position')
ylabel('Y Position')
zlabel('Vel. Cor. Signal')

Pos_X = Pos_X(:);
Pos_Y = Pos_Y(:);
Vx_Demo = Vx_Demo(:);
if FilterFundamentalSelect==0
    FilterFundamental=0;
end

if FilterFundamental==1
    SigX_Demo_TMP = bandstop(SigX_Demo(HighVelocityIndex),[fDriveX*.9,fDriveX*1.1],Fs);
    SigX_DemoDisp = SigX_Demo_TMP(:);
else
    SigX_DemoDisp = SigX_Demo(HighVelocityIndex);
    SigX_DemoDisp = SigX_DemoDisp(:);
end
HighVelocityIndex = abs(Vx_Demo)>(.3*max(Vx_Demo(:))); %This is just a vector with ones corresponding to the times where the velocity is above 0.3* the max
F = scatteredInterpolant(Pos_X(HighVelocityIndex),Pos_Y(HighVelocityIndex),SigX_DemoDisp./Vx_Demo(HighVelocityIndex)); %Vel COr

Interp_Data = F(X,Y);
ImageOut = Interp_Data;
% [Interp_Data_F0Cor] = SetZeroBC(Interp_Data);
figure,surf(Interp_Data);
xlabel( 'X Position')
ylabel('Y Position')
title('Surface plot of interpolated data')
% figure,surf(Interp_Data_F0Cor);
figure,imagesc(flipud(Interp_Data))
% figure,imagesc(flipud(Interp_Data_F0Cor))

axis image
xlabel( 'X Position')
ylabel('Y Position')
title('Image of interpolated data')

%%

figure
axes('Position',[.2 .1 .6 .55]);
plot(Pos_X(:),Pos_Y(:),'LineWidth',1.5)
title('FFP location over time','FontSize',13)
xlabel('X Position','FontSize',12)
ylabel('Y Position','FontSize',12)
axes('Position',[.15 .76 .7 .16]);

plot(tSig,Vx_Demo(:),'LineWidth',1.5)
title('FFP velocity over time (X only)','FontSize',13)
ylabel('FFP Vel.','FontSize',12)
xlabel('Time','FontSize',12)



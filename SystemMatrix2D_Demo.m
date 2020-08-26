function [ImageOut] = SystemMatrix2D_Demo(ImageIn,PlottingOn,SaveVid,varargin)
%% Introduction

% This function is a demonstration of the working principles behind
% system-maxtrix based Magnetic Particle Imaging (MPI). It makes a number of assumptions to
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
%   with path of "ImageIn__LissajousVid" if the function is given an image
%   path in, or default to "[Date]_LissajousVid"
NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments in (past "SaveVid") must be even')
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
ImageIn = ImageIn>max(max(ImageIn))/2; % If you want to make the object binary uncomment this




fDriveX = 25.5e3;%Hz
fDriveY = 24.e3;%Hz
Fs= 2.5e6;

x = -.02:.001:.02;%FOV meters
% x = linspace(-.04,.04,10);
y = x;
[X,Y] = meshgrid(x,y);
Grad = 5; %Gradient in Tesla per meter
BfflFunc = @(X) Grad.*X;
BfflX = BfflFunc(X); %B field of the FFL
BfflY = BfflFunc(Y); %B field of the FFL

t = 0:1/Fs:10e-3-(1/Fs);
BDriveX(1,1,:) = cos(2*pi*fDriveX*t);
BDriveX = repmat(BDriveX,length(x),length(x));
BDriveY(1,1,:) = cos(2*pi*fDriveY*t);
BDriveY = repmat(BDriveY,length(x),length(x));


L = length(t);
dt = 1/Fs;

f = Fs*(0:(L/2))/L;
%


NumHarmonics = round(f(end)/fDriveX);
Harmonics_Vec = [[(1:NumHarmonics)*fDriveX]'; [(1:NumHarmonics)*fDriveY]'];

Delta_f_Drive = abs(fDriveX-fDriveY);
Num_Mixing_Freqs = 10;
MixingFreqs_Vec = (-1*Num_Mixing_Freqs:Num_Mixing_Freqs)*Delta_f_Drive;
MixingFreqs_Vec_TMP = repmat(MixingFreqs_Vec,length(Harmonics_Vec),1);
Harmonics_Vec_TMP = repmat(Harmonics_Vec,1,length(MixingFreqs_Vec));

Acquired_Freqs = Harmonics_Vec_TMP+MixingFreqs_Vec_TMP;
Acquired_Freqs = Acquired_Freqs(:);
Harmonic_Locs = find(sum(mod(f,Acquired_Freqs)==0,1)>0);
BfflX = repmat(BfflX,1,1,length(t));
BX = BfflX+(1.1*Grad*x(end)*BDriveX);
BfflY = repmat(BfflY,1,1,length(t));
BY = BfflY+(1.1*Grad*x(end)*BDriveY);

ConcentrationMap = zeros(length(x)^2,1);
BMag = BX.^2+BY.^2;
% BMag(BMag==0) = 1;
clear BDriveX BDriveY BfflX BfflY

%%
Row_Demo = 10;
Col_Demo = 5;


figure('Position',[100 100 700 700])
Ax2 = axes('Position',[0.05 0.825 0.9 0.15]);
Ax4 = axes('Position',[0.05 0.62 0.9 0.125]);
Ax1 = axes('Position',[0.25 0.025 0.5 0.5]);

% Ax3 = axes('Position',[0.7 0.025 0.2 0.2]);
% imagesc(ParticleLoc_Image)
%     set(gca,'XTick',[], 'YTick', [])
% title('True Object')
Mx_Demo = Langevin(BMag(Col_Demo,Row_Demo,:),1).*BX(Col_Demo,Row_Demo,:)./(BMag(Col_Demo,Row_Demo,:));
My_Demo = Langevin(BMag(Col_Demo,Row_Demo,:),1).*BY(Col_Demo,Row_Demo,:)./(BMag(Col_Demo,Row_Demo,:));
SigX_Demo = diff(Mx_Demo(:));
SigY_Demo = diff(My_Demo(:));

if PlottingOn==1
    EndPt = length(t)-1;
    PlotRate = 3000; %This line controls how quickly the plotting occurs by skipping frames. If this is = 1, every frame is recorded, if it is equal to 10, every 10th frame is recorded
    %%
    VidCount=1;
    for jj = 1:PlotRate:EndPt
        
        
        axes(Ax1)
        cla
        imagesc(BX(:,:,jj).^2+BY(:,:,jj).^2)
        hold on
        plot(Row_Demo,Col_Demo,'ko','LineWidth',3)
        
        plot(Row_Demo,Col_Demo,'r.','LineWidth',5)
        
        title('Magnitude of B Field')
        caxis([0 .05])
        
        set(gca,'XTick',[], 'YTick', [])
        colormap parula
        axes(Ax2)
        cla
        % plot(t(1:jj),squeeze(sqrt(Mx(1:jj).^2+Mx(1:jj).^2)))
        plot(t(1:jj)*1000,squeeze(Mx_Demo(1:jj)),'r')
        hold on
        plot(t(1:jj)*1000,squeeze(My_Demo(1:jj)),'b')
        legend('X Magnetization','Y Magnetization')
        xlim([0 t(1000)*1000])
        xlabel('Time (milliseconds)','FontSize',12,'FontWeight','bold')
        
        axes(Ax4)
        cla
        % plot(t(1:jj),squeeze(sqrt(Mx(1:jj).^2+Mx(1:jj).^2)))
        plot(t(1:jj)*1000,SigX_Demo(1:jj),'r')
        hold on
        plot(t(1:jj)*1000,SigY_Demo(1:jj),'b')
        legend('X Signal','Y Signal')
        xlim([0 t(EndPt)*1000])
        xlabel('Time (milliseconds)','FontSize',12,'FontWeight','bold')
        
        
        if SaveVid==1
            frame_1(VidCount) = getframe(gcf);
            
            VidCount = VidCount+1;
        end
        
        
    end
    if SaveVid==1
        
        VideoPath = [ImagePath,'_LissajousVid'];
        
        v = VideoWriter(VideoPath,'MPEG-4');
        open(v)
        writeVideo(v,frame_1)
        close(v)
    end
end
%%
Ax =zeros(length(Harmonic_Locs(:)),length(x)*length(y));

Ay =zeros(length(Harmonic_Locs(:)),length(x)*length(y));

Build_A_Fig = figure('Position',[100 100 1000 800]);
AMatrixAxes = axes('Position',[.1 .1 .5 .8]);
CurrentSampleLocAxes = axes('Position',[.7 .55 .25 .35]);

FTAxes = axes('Position',[.7 .1 .25 .35]);
PlotRate2 = 3;
  VidCount=1;
for i = 1:length(x)*length(y) %For each position
    
    ConcentrationMap(i) = 1;
    ConcMap_ReShape = reshape(ConcentrationMap,length(x),length(x));
    [Col,Row] = find(ConcMap_ReShape==1);
    Mx = Langevin(BMag(Col,Row,:),1).*BX(Col,Row,:)./(BMag(Col,Row,:));
    My = Langevin(BMag(Col,Row,:),1).*BY(Col,Row,:)./(BMag(Col,Row,:));
    

    Ax(:,i) = Mag_2_A(squeeze(Mx),Harmonic_Locs);
    Ay(:,i) = Mag_2_A(squeeze(My),Harmonic_Locs);
    
    if mod(i,PlotRate2)==0
    axes(AMatrixAxes)
    cla
    imagesc(abs(vertcat(Ax,Ay)))
    title('The System Matrix')
    colorbar
    
    axes(CurrentSampleLocAxes)
    cla
    imagesc(ConcMap_ReShape)
    title('Test Sample Location')
    
    axes(FTAxes)
    cla
    Mag_X =abs(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
    Phase_X =angle(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
    plot(Mag_X)
    hold on
    Mag_Y =abs(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
    Phase_Y =angle(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
    plot(Mag_Y)
    title('Magnitude of Frequency Components')
    xlabel('Frequency')
    ylabel('Magnitude')
    legend('X Signal','Y Signal')
    
    sgtitle('Constructing the "System Matrix"','FontSize',16,'FontWeight','Bold')
    pause(.01)
    end
            
        if SaveVid==1
            frame_2(VidCount) = getframe(gcf);
            
            VidCount = VidCount+1;
        end
        
    disp(num2str(i))
    ConcentrationMap(i) = 0;
end
    if SaveVid==1
        
        VideoPath = [ImagePath,'_BuildSystemMatrix'];
        
        v = VideoWriter(VideoPath,'MPEG-4');
        open(v)
        writeVideo(v,frame_2)
        close(v)
    end

A = vertcat(Ax,Ay);

sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3) %Obnoxious beeping to alert me it is done
pause(1.5)
sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3) %More obnoxious beeping



%%


ConcentrationMap_ImageIn = imresize(ImageIn,size(X));
ConcentrationMap_ImageIn = cast(ConcentrationMap_ImageIn,'double');

Mx_Test = Langevin(BMag,ConcentrationMap_ImageIn).*BX./(BMag);
My_Test = Langevin(BMag,ConcentrationMap_ImageIn).*BY./(BMag);


Mag_Integrated_Test_X = squeeze(sum(sum(Mx_Test,2),1));
signal_Test_X = diff(Mag_Integrated_Test_X)';
Mag_Integrated_Test_Y = squeeze(sum(sum(My_Test,2),1));
signal_Test_Y = diff(Mag_Integrated_Test_Y)';

F_Test_X  = fft(signal_Test_X)';
F_Test_Y  = fft(signal_Test_Y)';

F_Test = [F_Test_X(Harmonic_Locs) ; F_Test_Y(Harmonic_Locs)] ;
 
TestRecon = pcg(A'*A,A'*F_Test,1e-9,300);
OneD_Image = abs(TestRecon);
figure,imagesc(reshape(OneD_Image,length(y),length(y)))
end
%%
function [A] = Mag_2_A(Mag,Harmonic_Locs)

signal(:) = diff(Mag)';
F(:)  = fft(signal);

A(:,1) = F(Harmonic_Locs);

end
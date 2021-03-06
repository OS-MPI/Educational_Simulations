function [ImageOut] = SystemMatrix2D_Demo(ImageIn,PlottingOn,SaveVid,varargin)
%% Introduction

% This function is a demonstration of the working principles behind
% system-maxtrix based Magnetic Particle Imaging (MPI). It makes a number of assumptions to
% simplify the system to make it easier to read and therefore learn from.
% Namely, the magnetic fields are assumed to be ideal (drive is
% homogeneous, and the FFP field is an ideal linear gradient).
% This function was written as a teaching tool and you should feel free to
% use it as such.

% The required inputs are:
%   ImageIn - A path to the object you want to simulate, or the object itself. It should be
%   a 2D matrix of numeric values
%   PlottingOn - Logical (1=on, 0 =off) that either enables or disables the
%   plotting of image as it is created.
%   SaveVid - if it is 1, the function will save a video of the plotting
%   with path of "ImageIn__LissajousVid" if the function is given an image
%   path in, or default to "[Date]_LissajousVid"

% Optionally, one can specify the X Drive frequency (fDriveX), Y Drive
% Frequency, sampling rate, sample particle position, and the animation
% plotting rate (how frequently images are plotted)
NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments in (past "SaveVid") must be even')
end
varargin = reshape(varargin,2,NumVarArgsIn/2);

fDriveX = 25e3;%Hz
fDriveY = 21e3;%Hz
Fs= 2.5e6;
PointSampleLocs = [20 30];
PlotRate = 1; %This line controls how quickly the plotting occurs by skipping frames. If this is = 1, every frame is recorded, if it is equal to 10, every 10th frame is recorded

for i = 1:NumVarArgsIn/2
    if strcmp(varargin(1,i),'AddNoise')
        AddNoise = cell2mat(varargin(2,i));
        AddNoiseSelected = 1; %Indicate the user selected to add noise or not
    elseif exist('AddNoiseSelected','var')
    else
        AddNoiseSelected = 0;
    end

    
    if strcmp(varargin(1,i),'PointSampleLocs')
        PointSampleLocs = cell2mat(varargin(2,i));
        PointSampleSelect = 1; 
    end
    if strcmp(varargin(1,i),'fs')
        fs = cell2mat(varargin(2,i));
    end
    
    if strcmp(varargin(1,i),'fDriveX')
        fDriveX = cell2mat(varargin(2,i));
    end
    
    if strcmp(varargin(1,i),'fDriveY')
        fDriveY = cell2mat(varargin(2,i));
    end
    
    if strcmp(varargin(1,i),'PlotRate')
        PlotRate = cell2mat(varargin(2,i));
    end
    
    if strcmp(varargin(1,i),'SoundOn')
        SoundOn = cell2mat(varargin(2,i)); %Any Value of SoundOn besides 0 will cause beeping when completed
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






PointsPerDrivePrd = round(Fs/fDriveX);

x = -.02:.001:.02;%FOV meters
% x = linspace(-.04,.04,10);
y = x;
[X,Y] = meshgrid(x,y);
Grad = 5; %Gradient in Tesla per meter
BfflFunc = @(X) Grad.*X;
BfflX = BfflFunc(X); %B field of the FFL
BfflY = BfflFunc(Y); %B field of the FFL

t = 0:1/Fs:30e-3-(1/Fs);
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

ConcentrationVector = zeros(length(x)*length(y),1);
BMag = sqrt(BX.^2+BY.^2);
% BMag(BMag==0) = 1;
clear BDriveX BDriveY BfflX BfflY
ConcentrationMap = zeros(size(BMag(:,:,1)));


if PointSampleLocs(1)>size(ConcentrationMap,2)
    PointSampleLocs(1) = size(ConcentrationMap,2);
end
if PointSampleLocs(2)>size(ConcentrationMap,1)
    PointSampleLocs(2) = size(ConcentrationMap,1);
end

PointSampleLocs(2) = size(ConcentrationMap,1)-PointSampleLocs(2)+1;

ConcentrationMap(PointSampleLocs(2),PointSampleLocs(1)) = 1;
tmpCoords = [PointSampleLocs(2),PointSampleLocs(1)];


%%



figure('Position',[100 100 700 800])
AnimationFromBox = annotation('textbox',[.0,0,.2,.025],'String','Animation from OS-MPI.GitHub.io','FitBoxToText','on');
sgtitle('Lissajous excitation MPI Signal generation')
Ax2 = axes('Position',[0.05 0.8 0.9 0.1]);
Ax4 = axes('Position',[0.05 0.62 0.9 0.1]);
Ax1 = axes('Position',[0.25 0.1 0.5 0.4]);

% Ax3 = axes('Position',[0.7 0.025 0.2 0.2]);
% imagesc(ParticleLoc_Image)
%     set(gca,'XTick',[], 'YTick', [])
% title('True Object')
Mx_Demo = Langevin(BMag(PointSampleLocs(2),PointSampleLocs(1),:),1).*BX(PointSampleLocs(2),PointSampleLocs(1),:)./(BMag(PointSampleLocs(2),PointSampleLocs(1),:));
My_Demo = Langevin(BMag(PointSampleLocs(2),PointSampleLocs(1),:),1).*BY(PointSampleLocs(2),PointSampleLocs(1),:)./(BMag(PointSampleLocs(2),PointSampleLocs(1),:));
TotalMag = sqrt(squeeze(Mx_Demo(:)).^2+squeeze(My_Demo(:)).^2);
% MaxMag = max([max(Mx_Demo),max(My_Demo)]);
MaxMag = max(TotalMag(:));
MinMag = min([min(Mx_Demo),min(My_Demo)]);




SigX_Demo = diff(Mx_Demo(:));
SigY_Demo = diff(My_Demo(:));
MaxSig = max([max(SigX_Demo),max(SigY_Demo)]);
MinSig = min([min(SigX_Demo),min(SigY_Demo)]);



if PlottingOn==1
    EndPt = length(t)-1;
    
    TimeVecDrivePrds = 10;
    TimeVecNumPoints = TimeVecDrivePrds*PointsPerDrivePrd;
    
    %%
    VidCount=1;
    for jj = 1:PlotRate:EndPt
        
        
        axes(Ax1)
        cla
        imagesc(x,y,BX(:,:,jj).^2+BY(:,:,jj).^2)
        %         colorbar
        hold on
        plot(x(PointSampleLocs(1)),y(PointSampleLocs(2)),'ko','LineWidth',3)
        
        
        quiver(x(PointSampleLocs(1)),y(PointSampleLocs(2)),Mx_Demo(jj),My_Demo(jj),5e-3,'r','LineWidth',3,'MaxHeadSize',10)
        
        
        xlabel('X Position')
        ylabel('Y Position')
        title('Magnitude of B Field')
        caxis([0 .05])
        legend('Nanoparticles','Magnetization Vector')
        set(gca,'YDir','normal')
        %         set(gca,'XTick',[], 'YTick', [])
        colormap parula
        axes(Ax2)
        cla
        %         plot(t(1:jj),squeeze(sqrt(Mx(1:jj).^2+Mx(1:jj).^2)))
        plot(t(1:jj)*1000,squeeze(Mx_Demo(1:jj)),'r')
        hold on
        plot(t(1:jj)*1000,squeeze(My_Demo(1:jj)),'b')
        plot(t(1:jj)*1000,TotalMag(1:jj),'g')
        if jj<=ceil(TimeVecNumPoints/2)+1
            xlim([0 t(TimeVecNumPoints)*1000])
        elseif jj>ceil(TimeVecNumPoints/2) && jj<EndPt-floor(TimeVecNumPoints/2)
            xlim([t(jj-floor(TimeVecNumPoints/2)) t(jj+floor(TimeVecNumPoints/2))]*1000)
        else
            xlim([t(end-floor(TimeVecNumPoints)) t(end)]*1000)
        end
        legend('X Magnetization','Y Magnetization','Net Magnetization')
        ylim([MinMag MaxMag]*1.1)
        xlabel('Time (milliseconds)','FontSize',12,'FontWeight','bold')
        
        axes(Ax4)
        cla
        % plot(t(1:jj),squeeze(sqrt(Mx(1:jj).^2+Mx(1:jj).^2)))
        plot(t(1:jj)*1000,SigX_Demo(1:jj),'r')
        hold on
        plot(t(1:jj)*1000,SigY_Demo(1:jj),'b')
        legend('X Signal','Y Signal')
        ylim([MinSig MaxSig])
        set(Ax4,'XLim',Ax2.XLim)
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
AMatrixAxes = axes('Position',[.1 .25 .5 .5]);
CurrentSampleLocAxes = axes('Position',[.7 .55 .25 .35]);

FTAxes = axes('Position',[.7 .1 .25 .35]);
PlotRate2 = 3;
VidCount=1;
for i = 1:length(x)*length(y) %For each position
    
    ConcentrationVector(i) = 1;
    ConcMap_ReShape = reshape(ConcentrationVector,length(y),length(x));
    [Col,Row] = find(ConcMap_ReShape==1);
    Mx = Langevin(BMag(Col,Row,:),1).*BX(Col,Row,:)./(BMag(Col,Row,:));
    My = Langevin(BMag(Col,Row,:),1).*BY(Col,Row,:)./(BMag(Col,Row,:));
    
    
    Ax(:,i) = Mag_2_A(squeeze(Mx),Harmonic_Locs);
    Ay(:,i) = Mag_2_A(squeeze(My),Harmonic_Locs);
    
    if mod(i,PlotRate2)==0
        axes(AMatrixAxes)
        cla
        imagesc(abs(Ax))
        title({'The System Matrix';'Single Rx Coil'})
        ylabel('Frequency Component')
        xlabel('Particle Index')
        colorbar
        
        axes(CurrentSampleLocAxes)
        cla
        imagesc(x,y,ConcMap_ReShape)
        set(gca,'YDir','normal')
        xlabel('X Position')
        ylabel('Y Position')
        title('Test Sample Location')
        
        axes(FTAxes)
        cla
        Mag_X =abs(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
        Phase_X =angle(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
        plot(Mag_X)
        hold on
        Mag_Y =abs(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
        Phase_Y =angle(Ax(1:round(length(Harmonic_Locs)/2)-1,i));
        %     plot(Mag_Y)
        title('Magnitude of Frequency Components')
        xlabel('Frequency')
        ylabel('Magnitude')
        %     legend('X Signal')
        
        sgtitle('Constructing the "System Matrix"','FontSize',16,'FontWeight','Bold')
        pause(.01)
    end
    
    if SaveVid==1
        frame_2(VidCount) = getframe(gcf);
        
        VidCount = VidCount+1;
    end
    
    disp(num2str(i))
    ConcentrationVector(i) = 0;
end
if SaveVid==1
    
    VideoPath = [ImagePath,'_BuildSystemMatrix'];
    
    v = VideoWriter(VideoPath,'MPEG-4');
    open(v)
    writeVideo(v,frame_2)
    close(v)
end

A = vertcat(Ax,Ay);

% sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3) %Obnoxious beeping to alert me it is done
% pause(1.5)
% sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3) %More obnoxious beeping

if exist('SoundOn','var')
    if SoundOn~=0
        MorseBeep('.... ..');
    end
end

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
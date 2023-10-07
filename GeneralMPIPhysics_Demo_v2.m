function [] = GeneralMPIPhysics_Demo_v2(SaveVid,varargin)
%% Function Inputs/ help

% SaveVid is a binary option, if you want to save a video (1) or not (0)

% For Varargin:
% - FileName: is the file name that will be saved
% - VidFormat: either "gif" or "mp4"
% - MoveSample: Binary 1/0 if you want the sample to move
% - SoundOn: Binary 1/0 if you want it to make a beep when it finishes

%%

addpath([pwd,'\Called Functions'])

VidFormat='mp4'; %Default. WIll change if user enters something else
filename = [date(),'MPI_Physics_Vid.mp4']; %Default value. Will change if user enters something else when starting function
MoveSample = 0; %Default
yloc = 0;
yloc_FFP =0;
gifloops = inf;
NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments inputted must be even')
elseif NumVarArgsIn==0
    
    filename = [date(),'MPI_Physics_Vid_v2.mp4'];
    
end

varargin = reshape(varargin,2,NumVarArgsIn/2);

for i = 1:NumVarArgsIn/2
    if strcmpi(varargin(1,i),'FileName')
        filename = cell2mat(varargin(2,i));
    end
    
    if strcmpi(varargin(1,i),'VidFormat')
        VidFormat = cell2mat(varargin(2,i));
        
        VidFormat(VidFormat=='.')=[]; %Making sure there is no period in the format. So ".GIF" goes to "GIF", Just for consistency
        
        if not(strcmpi(VidFormat,'mp4')|strcmpi(VidFormat,'gif'))
            error('Video format must be mp4 or gif')
        end 
        filename = [date(),'MPI_Physics_Vid.',VidFormat];

    end
    
    if strcmpi(varargin(1,i),'MoveSample') %Selecting MoveSample as 1(Yes) or 0(No) moves the sample along the axis to show how the signal changes with biasing
        MoveSample = cell2mat(varargin(2,i));
    else
        
        
    end
    if strcmpi(varargin(1,i),'SoundOn')
        SoundOn = cell2mat(varargin(2,i));
    end
    if strcmpi(varargin(1,i),'MorseIn')
        MorseIn = cell2mat(varargin(2,i));
    end
    if strcmpi(varargin(1,i),'YLoc')
        yloc = cell2mat(varargin(2,i));
    end
    if strcmpi(varargin(1,i),'YLoc_FFP')
        yloc_FFP = cell2mat(varargin(2,i));
    end
    
    if strcmpi(varargin(1,i),'GIFLoops')
        gifloops = cell2mat(varargin(2,i));
    end


    
end

filename_HasExt = sum(filename=='.')>0;
if filename_HasExt==1
    ExtensionStart = find(filename=='.');
    VidFormat = filename(ExtensionStart+1:end);
else
    filename = [filename,'.',VidFormat];
end

SavePath = uigetdir();
filename = [SavePath,'\',filename];
isGIF = strcmpi(VidFormat,'gif');
isMP4 = strcmpi(VidFormat,'mp4');


fDrive = 25e3;%Hz
Fs= 2e6; %Hz If this is too low, you start to get discretization artifacts in the figures. Keep above ~50x fDrive

x = -.025:.001:.025;%FOV meters
y = x;
y = y-yloc_FFP;
x_downSamp = x(1:5:end);
y_downSamp = x_downSamp;
[X,Y] = meshgrid(x_downSamp,y_downSamp);

XZero_Ind = find(x==0);


if exist('MoveSample','var')
    if MoveSample==1
        X_EndIndex=length(x);
    else
        X_EndIndex = XZero_Ind;
    end
else %The default if you do not select motion or no is to not move
    X_EndIndex = XZero_Ind;
end

Grad = 2; %Gradient in Tesla per meter
BfflFunc = @(x) Grad.*x;
Bffl = BfflFunc(x); %B field of the FFL
B_FFP_X= BfflFunc(X); %B field of the FFP
B_FFP_Y = BfflFunc(Y); %B field of the FFP

t = 0:1/Fs:3e-3-(1/Fs);
t_Up = 0:1/(Fs*4):3e-3-(1/Fs);
dt = 1/Fs;


Count = 0;
GIF_Count = 1;
% for XPos_Ind = XZero_Ind:1:length(x)
for XPos_Ind = XZero_Ind:X_EndIndex
    Count = Count +1;
    
    XPos = x(XPos_Ind);
    
    B = Grad*x(end)/1.5*cos(2*pi*fDrive*t);
    B_Big  = 2*Grad*x(end)*cos(2*pi*fDrive*t_Up);
    
    B_Bias = sqrt((B+ Grad*XPos).^2 + (Grad*(yloc-yloc_FFP)).^2);
    if (yloc-yloc_FFP)==0
        B_Bias = B_Bias.* sign(B);
    end
    M = Langevin(B,1);
    M_Big = Langevin(B_Big,1);
    M_Bias = Langevin(B_Bias,1);
    
    Signal = diff(M);
    FT_Sig = fft(Signal(1:end-1));
    
    L_FT = length(FT_Sig);
    P1_FT = abs(FT_Sig(1:L_FT/2+1));
    f = Fs*(0:L_FT/2)/L_FT;
    %%
    
    if Count==1
        h= figure('Position',[100 100 1300 850]);set(gcf, 'Color', 'White');
        CitationBox = annotation('textbox',[.0,0,.1,.025],'String',{'Adapted from Gleich and Weizenecker, 2005 Fig. 1'},'FitBoxToText','on');
        AnimationFromBox = annotation('textbox',[.0,0,.2,.025],'String','Animation from OS-MPI.GitHub.io','FitBoxToText','off');
        drawnow
        AnimationFromBox.Position(1) = 1-AnimationFromBox.Position(3);
        AnimationFromBox.Position(2) = CitationBox.Position(2);
        AnimationFromBox.Position(4) = CitationBox.Position(4);
        P1 = subplot(2,3,1);
        P1.Position(2) = P1.Position(2);
        P1.Position(4) = P1.Position(4);
        
        P2 = subplot(2,3,2);
        P2.Position(2) = P2.Position(2);
        P2.Position(4) = P2.Position(4);
        
        P3 = subplot(2,3,3);
        P3.Position(2) = P3.Position(2);
        P3.Position(4) = P3.Position(4);
        
        P4 = subplot(2,3,4);
        P4.Position(2) = P4.Position(2);
        P4.Position(4) = P4.Position(4);
        
        P5 = subplot(2,3,5);
        P5.Position(2) = P5.Position(2);
        P5.Position(4) = P5.Position(4);
        
        
        P6 = subplot(2,3,6);
        P6.Position(2) = P6.Position(2);
        P6.Position(4) = P6.Position(4);
        
%         TMP_Ax = axes('Position',[0 0 1 .25]);
%         rectangle('Position',[0 0 1 .125],'FaceColor',[.85 .85 .85])
%         TMP_Ax.Color = 'None';
%         set(gca,'XTick',[], 'YTick', [])
%         TMP_Ax.YColor = 'None';
%         TMP_Ax.XColor = 'None';
        
        
        AxGap14 = axes('Position',[P1.Position(1) P4.Position(2)+P4.Position(4) P4.Position(3) P1.Position(2)-(P4.Position(2)+P4.Position(4))]);
        Ax14Over = axes('Position',P1.Position);
        
        AxGap12 = axes('Position',[P1.Position(1)+P1.Position(3) P1.Position(2) P2.Position(1)-(P1.Position(1)+P1.Position(3)) P1.Position(4)]);
        Ax12Over = axes('Position',P1.Position);
        
%         FFPAx = axes('Position',[P1.Position(1) .05 (P6.Position(1)+P6.Position(3)-P1.Position(1)) .125]);
        
    end
    
    BPick = B_Bias;
    MPick = M_Bias;
    Signal = diff(MPick);
    
    FT_Sig = fft(Signal(1:end-1));
    
    L_FT = length(FT_Sig);
    P1_FT = abs(FT_Sig(1:L_FT/2+1));
    if Count==1
        P1_FT_Orig = P1_FT;
    end
    f = Fs*(0:L_FT/2)/L_FT;
    if mod(Count,2)==0
        StartI = round(Fs/fDrive/2)+1;
        EndI = round(Fs/fDrive/2)*2;
    else
        StartI = 1;
        EndI = round(Fs/fDrive/2);
        if MoveSample==0
            EndI = round(Fs/fDrive);
        end
    end
    for ii = StartI:EndI
        
        if ii==StartI && Count==1
            axes(P1)
            PlotMBLine = plot(B_Big(1:400),M_Big(1:400),'k','LineWidth',3);
            hold on
            Dot = plot(BPick(ii),MPick(ii),'bo','MarkerSize',10,'LineWidth',2);
            hold off
            xlim([min(B_Big) max(B_Big)])
            xlabel('External Magnetic Field','FontSize',14,'FontWeight','bold')
            title('Particle Magnetization Curve','FontSize',14,'FontWeight','bold')
            ylabel({'SPION';'Magnetization'},'FontSize',14,'FontWeight','bold')
        else
            set(PlotMBLine,'XData',B_Big(1:400))
            set(PlotMBLine,'YData',M_Big(1:400))
            set(Dot,'XData',BPick(ii))
            set(Dot,'YData',MPick(ii))
            
        end
        
        
        
        
        
        
        
        if ii==StartI && Count==1
            axes(P4)
            FieldOverTime = plot(BPick(1+ii:200+ii),t(1+ii:200+ii),'Color',[30,119,20]/255,'LineWidth',3);
            ylim([ t(1+ii) t(200+ii)])
            xlim([min(B_Big) max(B_Big)])
            axis ij
            ylabel('Time (seconds)','FontSize',14,'FontWeight','bold')
            xlabel('External Magnetic Field','FontSize',14,'FontWeight','bold')
        else
            set(FieldOverTime,'XData',BPick(1+ii:200+ii))
            set(FieldOverTime,'YData',t(1+ii:200+ii))
            
            set(P4,'ylim',([ t(1+ii) t(200+ii)]))
            set(P4,'xlim',[min(B_Big) max(B_Big)])
            
            
        end
        
        if ii==StartI&& Count==1
            axes(P2)
            MagOverTime = plot(t(1+ii:200+ii),MPick(1+ii:200+ii),'k','LineWidth',3);
            xlim([t(1+ii) t(200+ii)])
            ylim(P1.YLim)
            xlabel('Time (seconds)','FontSize',14,'FontWeight','bold')
            title('Magnetization over time','FontSize',14,'FontWeight','bold')
        else
            set(MagOverTime,'YData',MPick(1+ii:200+ii))
            set(MagOverTime,'XData',t(1+ii:200+ii))
            
            set(P2,'xlim',([ t(1+ii) t(200+ii)]))
            
        end
        
        if ii==StartI&& Count==1
            axes(P3)
            SigOverTime = plot(t(1+ii:200+ii),Signal(1+ii:200+ii),'Color',[33,204,192]/255,'LineWidth',3);
            xlim([t(1+ii) t(200+ii)])
            ylim([-0.5 0.5])
            xlabel('Time (seconds)','FontSize',14,'FontWeight','bold')
            ylabel('dM/dt','FontSize',14,'FontWeight','bold')
            title('Signal over time','FontSize',14,'FontWeight','bold')
        else
            set(SigOverTime,'YData',Signal(1+ii:200+ii))
            set(SigOverTime,'XData',t(1+ii:200+ii))
            
            set(P3,'xlim',([ t(1+ii) t(200+ii)]))
            
        end
        
        if ii==StartI&& Count==1
            axes(P5)
            
            FieldQuiver = quiver(X,Y,-B_FFP_X+B(1+ii),-B_FFP_Y+Grad*yloc_FFP,1,'Color',[30,119,20]/255,'LineWidth',2);
            
            ylim([min(x) max(x)])
            xlim([min(x) max(x)])
            axis image
            
            title('External field','FontSize',14,'FontWeight','bold')
            hold on
            ParticlePlot0 = plot(x(XPos_Ind),yloc,'ko','LineWidth',4,'MarkerSize',1);
            ParticlePlot1 = plot(x(XPos_Ind),yloc,'ko','LineWidth',4,'MarkerSize',10);
            ParticlePlot2 = plot(x(XPos_Ind),yloc,'ro','LineWidth',2);
            FFPPlot = plot(B(1+ii)/Grad,yloc_FFP,'mo','LineWidth',2,'MarkerSize',30);
            
            SPIONField = quiver(x(XPos_Ind),yloc,B(1+ii),-BfflFunc(yloc-yloc_FFP),0.1,'r','LineWidth',3);
%         set(gca,'XTick',[1:5:length(x)], 'YTick', [])
%         FFPAx.XTickLabel = {x(1:5:length(x))};
%         title('B-Field along 1D Line')
%         xlabel('Position (meters)')
%         Lgd = legend('','Particle location');
%         Lgd.Location = 'northeastoutside';
        else
            set(FieldQuiver,'UData',-B_FFP_X+B(1+ii))
            set(SPIONField,'UData',B(1+ii))
            set(SPIONField,'XData',x(XPos_Ind))
%             set(FieldQuiver,'YData',t(1+ii:200+ii))
            set(ParticlePlot0,'XData',x(XPos_Ind))
            set(ParticlePlot1,'XData',x(XPos_Ind))
            set(ParticlePlot2,'XData',x(XPos_Ind))
            set(FFPPlot,'XData',B(1+ii)/Grad)
            set(P5,'ylim',[min(x) max(x)])
            set(P5,'xlim',[min(x) max(x)])
            
            
        end
        
        if ii==StartI&& Count==1
            axes(P6)
            FTPlot = plot(f/fDrive,P1_FT,'Color',[33,204,192]/255,'LineWidth',3);
            ylim([0 max(P1_FT_Orig)])
            xlabel('Frequencies (f/f_{Drive})','FontSize',14,'FontWeight','bold')
            ylabel('Magnitude','FontSize',14,'FontWeight','bold')
        else
            set(FTPlot,'YData',P1_FT)
        end
        
        
        
%         axes(FFPAx)
%         imagesc(repmat(abs(Bffl+B(ii)),3,1))
%         caxis([0 max(abs(B))*1.25])
%         if ii==1 && Count==1
%             CBR = colorbar;
%             CBR.Label.String = {'External Field'; 'Magnitude'};
%             CBR.Label.FontWeight = 'Bold';
%             CBR.Location = 'westoutside';
%         end
%         hold on
%         plot(XPos_Ind,2,'ko','LineWidth',4)
%         %     plot(XPos_Ind,2,'bo','LineWidth',1)
%         set(gca,'XTick',[1:5:length(x)], 'YTick', [])
%         FFPAx.XTickLabel = {x(1:5:length(x))};
%         title('B-Field along 1D Line')
%         xlabel('Position (meters)')
%         Lgd = legend('Particle location');
%         Lgd.Location = 'northeastoutside';
        
        
        
        
        if ii==StartI&& Count==1
            
            axes(AxGap14)
            Axes14GapPlot = plot([BPick(ii) BPick(ii)],[0 1],'b--');
            AxGap14.Color = 'None';
            set(gca,'XTick',[], 'YTick', [])
            AxGap14.YColor = 'None';
            AxGap14.XColor = 'None';
            xlim(P1.XLim)
            ylim([0 1])
            
        else
            set(Axes14GapPlot,'XData',[BPick(ii) BPick(ii)])
        end
        
        if ii==StartI && Count==1
            axes(Ax14Over)
            Axes14OverFigPlot = plot([BPick(ii) BPick(ii)],[P1.YLim(1) MPick(ii) ],'b--');
            Ax14Over.Color = 'None';
            set(gca,'XTick',[], 'YTick', [])
            Ax14Over.YColor = 'None';
            Ax14Over.XColor = 'None';
            xlim(P1.XLim)
            ylim(P1.YLim)
        else
            set(Axes14OverFigPlot,'XData',[BPick(ii) BPick(ii)])
            set(Axes14OverFigPlot,'YData',[P1.YLim(1) MPick(ii) ])
        end
        
        
        if ii==StartI && Count==1
            axes(AxGap12)
            Axes12GapPlot = plot([0 1],[MPick(ii) MPick(ii)],'b--');
            AxGap12.Color = 'None';
            set(gca,'XTick',[], 'YTick', [])
            AxGap12.YColor = 'None';
            AxGap12.XColor = 'None';
            ylim(P1.YLim)
            xlim([0 1])
        else
            set(Axes12GapPlot,'YData',[MPick(ii) MPick(ii)])
            
        end
        
        
        
        if ii==StartI && Count==1
            axes(Ax12Over)
            
            Axes12OverPlot = plot([BPick(ii) P1.XLim(2) ],[MPick(ii) MPick(ii) ],'b--');
            Ax12Over.Color = 'None';
            set(gca,'XTick',[], 'YTick', [])
            Ax12Over.YColor = 'None';
            Ax12Over.XColor = 'None';
            xlim(P1.XLim)
            ylim(P1.YLim)
        else
            set(Axes12OverPlot,'XData',[BPick(ii) P1.XLim(2) ])
            set(Axes12OverPlot,'YData',[MPick(ii) MPick(ii) ])
        end
        
        
        drawnow
        %         pause(.001)
        if SaveVid==1
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
        end
        if isGIF && SaveVid==1
            if GIF_Count == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',gifloops,'DelayTime',1/32);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/32);
            end
        end
        disp([num2str(ii),',' num2str(GIF_Count)])
        
        frame_1(GIF_Count) = getframe(gcf);
        clear tmp
        GIF_Count = GIF_Count+1;
    end
    
%     sound(sawtooth(2*pi*.3e3*(0:1/10e3:1)).*sin(2*pi*.3e3*(0:1/10e3:1)),10e3) %Obnoxious beeping to alert me it is done
end
if isMP4 && SaveVid==1
    
    v = VideoWriter(filename,'MPEG-4');
    open(v)
    writeVideo(v,frame_1)
    close(v)
end

if exist('SoundOn','var')
    if SoundOn==1
    if not(exist('MorseIn','var'))
        MorseBeep('--  .-.-  ..')
    else
        if MorseIn==1
            MorseBeep(MorseIn)
        else
            MorseBeep('--  .-.-  ..')
        end
    end
    end
end


end

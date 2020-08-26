clear all
close all
% filename = 'GleichFig_Test_4.mp4';
filename = 'GleichFig_Test_4Center_TMP.gif';

fDrive = 25e3;%Hz
Fs= 2e6;

x = -.025:.001:.025;%FOV meters
XZero_Ind = find(x==0);
Grad = 2; %Gradient in Tesla per meter
BfflFunc = @(x) Grad.*x;
Bffl = BfflFunc(x); %B field of the FFL
t = 0:1/Fs:3e-3-(1/Fs);
t_Up = 0:1/(Fs*4):3e-3-(1/Fs);
dt = 1/Fs;



C = zeros(size(x));
C(x==0) = 1;
signal = zeros(length(t)-1,length(Bffl));
Count = 0;
GIF_Count = 1;
% for XPos_Ind = XZero_Ind:1:length(x)
for XPos_Ind = XZero_Ind
    Count = Count +1;
    F = zeros(length(t)-1,length(x));
    % XPos_Ind = length(x)-5;
    % XPos_Ind = XZero_In;
    XPos = x(XPos_Ind);
    
    B = Grad*x(end)/1.5*cos(2*pi*fDrive*t);
    B_Big  = 2*Grad*x(end)*cos(2*pi*fDrive*t_Up);
    
    B_Bias = B+ Grad*XPos;
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
        h= figure('Position',[100 100 1300 850]);
        P1 = subplot(2,3,1);
        P1.Position(2) = P1.Position(2)+.1;
        P1.Position(4) = P1.Position(4)-.1;
        
        P2 = subplot(2,3,2);
        P2.Position(2) = P2.Position(2)+.1;
        P2.Position(4) = P2.Position(4)-.1;
        
        P3 = subplot(2,3,3);
        P3.Position(2) = P3.Position(2)+.1;
        P3.Position(4) = P3.Position(4)-.1;
        
        P4 = subplot(2,3,4);
        P4.Position(2) = P4.Position(2)+.2;
        P4.Position(4) = P4.Position(4)-.1;
        
        
        P6 = subplot(2,3,6);
        P6.Position(2) = P6.Position(2)+.2;
        P6.Position(4) = P6.Position(4)-.1;
        
        TMP_Ax = axes('Position',[0 0 1 .25]);
        rectangle('Position',[0 0 1 .125],'FaceColor',[.85 .85 .85])
        TMP_Ax.Color = 'None';
        set(gca,'XTick',[], 'YTick', [])
        TMP_Ax.YColor = 'None';
        TMP_Ax.XColor = 'None';
        
        
        AxGap14 = axes('Position',[P1.Position(1) P4.Position(2)+P4.Position(4) P4.Position(3) P1.Position(2)-(P4.Position(2)+P4.Position(4))]);
        Ax14Over = axes('Position',P1.Position);
        
        AxGap12 = axes('Position',[P1.Position(1)+P1.Position(3) P1.Position(2) P2.Position(1)-(P1.Position(1)+P1.Position(3)) P1.Position(4)]);
        Ax12Over = axes('Position',P1.Position);
        
        FFPAx = axes('Position',[P1.Position(1) .05 (P6.Position(1)+P6.Position(3)-P1.Position(1)) .125]);
        
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
        StartI = 41;
        EndI = 80;
    else
        StartI = 1;
        EndI = 80;
    end
    for ii = StartI:EndI
        
        axes(P1)
        plot(B_Big(1:400),M_Big(1:400))
        hold on
        Dot = plot(BPick(ii),MPick(ii),'bo');
        hold off
        xlim([min(B_Big) max(B_Big)])
        xlabel('External Magnetic Field')
        title('Particle Magnetization Curve')
        ylabel('Magnetization')
        
        axes(P4)
        plot(BPick(1+ii:200+ii),t(1+ii:200+ii))
        ylim([ t(1+ii) t(200+ii)])
        xlim([min(B_Big) max(B_Big)])
        axis ij
        ylabel('Time (seconds)')
        xlabel('External Magnetic Field')
        
        axes(P2)
        plot(t(1+ii:200+ii),MPick(1+ii:200+ii))
        xlim([t(1+ii) t(200+ii)])
        ylim(P1.YLim)
        xlabel('Time (seconds)')
        title('Magnetization over time')
        
        axes(P3)
        plot(t(1+ii:200+ii),Signal(1+ii:200+ii))
        xlim([t(1+ii) t(200+ii)])
        ylim([-0.5 0.5])
        xlabel('Time (seconds)')
        ylabel('dM/dt')
        title('Signal over time')
        
        axes(P6)
        plot(f/fDrive,P1_FT)
        ylim([0 max(P1_FT_Orig)])
        xlabel('Frequencies (f/f_{Drive})')
        ylabel('Magnitude')
        
        
        axes(FFPAx)
        imagesc(repmat(abs(Bffl+B(ii)),3,1))
        caxis([0 max(abs(B))*1.25])
        if ii==1 && Count==1
            CBR = colorbar;
            CBR.Label.String = {'External Field'; 'Magnitude'};
            CBR.Label.FontWeight = 'Bold';
            CBR.Location = 'westoutside';
        end
        hold on
        plot(XPos_Ind,2,'ko','LineWidth',4)
        %     plot(XPos_Ind,2,'bo','LineWidth',1)
        set(gca,'XTick',[1:5:length(x)], 'YTick', [])
        FFPAx.XTickLabel = {x(1:5:length(x))};
        title('B-Field along 1D Line')
        xlabel('Position (meters)')
        Lgd = legend('Particle location');
        Lgd.Location = 'northeastoutside';
        
        
        
        axes(AxGap14)
        
        
        
        tmp = plot([BPick(ii) BPick(ii)],[0 1],'b--');
        AxGap14.Color = 'None';
        set(gca,'XTick',[], 'YTick', [])
        AxGap14.YColor = 'None';
        AxGap14.XColor = 'None';
        xlim(P1.XLim)
        ylim([0 1])
        
        
        axes(Ax14Over)
        
        
        plot([BPick(ii) BPick(ii)],[P1.YLim(1) MPick(ii) ],'b--')
        Ax14Over.Color = 'None';
        set(gca,'XTick',[], 'YTick', [])
        Ax14Over.YColor = 'None';
        Ax14Over.XColor = 'None';
        xlim(P1.XLim)
        ylim(P1.YLim)
        
        
        
        axes(AxGap12)
        tmp2 = plot([0 1],[MPick(ii) MPick(ii)],'b--');
        AxGap12.Color = 'None';
        set(gca,'XTick',[], 'YTick', [])
        AxGap12.YColor = 'None';
        AxGap12.XColor = 'None';
        ylim(P1.YLim)
        xlim([0 1])
        
        
        axes(Ax12Over)
        
        plot([BPick(ii) P1.XLim(2) ],[MPick(ii) MPick(ii) ],'b--')
        Ax12Over.Color = 'None';
        set(gca,'XTick',[], 'YTick', [])
        Ax12Over.YColor = 'None';
        Ax12Over.XColor = 'None';
        xlim(P1.XLim)
        ylim(P1.YLim)
        
        
        
        pause(.001)
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if GIF_Count == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/32);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/32);
        end
        
        disp([num2str(ii),',' num2str(GIF_Count)])
      
        frame_1(GIF_Count) = getframe(gcf);
        clear tmp
          GIF_Count = GIF_Count+1;
    end
    
    %     clf
    sound(sawtooth(2*pi*.3e3*(0:1/10e3:1)).*sin(2*pi*.3e3*(0:1/10e3:1)),10e3) %Obnoxious beeping to alert me it is done
end
%  v = VideoWriter(filename,'MPEG-4');
%         open(v)
%         writeVideo(v,frame_1)
%         close(v)

sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3) %Obnoxious beeping to alert me it is done
pause(1.5)
sound(sawtooth(2*pi*1e3*(0:1/10e3:1)).*sin(2*pi*1e3*(0:1/10e3:1)),10e3)
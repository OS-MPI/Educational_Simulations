function [] = MorseBeep(StringIn,varargin)


NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments inputted must be even')
end
varargin = reshape(varargin,2,NumVarArgsIn/2);
fs = 20e3;% Hz sampling frequency of sound
SoundFreq = 440; %Hz, Frequency of sound played
DotTime = 0.1; %How long a dot is
for i = 1:NumVarArgsIn/2
    if strcmpi(varargin(1,i),'Fs')
        fs = cell2mat(varargin(2,i));
    end
    
    if strcmpi(varargin(1,i),'SoundFreq')
        SoundFreq = cell2mat(varargin(2,i));
    end
    if strcmpi(varargin(1,i),'DotTime')
        DotTime = cell2mat(varargin(2,i));
    end
    
end

t = 0:1/fs:(DotTime-1/fs);
t = t(:); %Making a column



if sum((StringIn=='.')+(StringIn=='-')+(StringIn==' '))<length(StringIn)
    error('Input must be either "." or "-" or " ")')
end
SoundVec = [];
for i=1:length(StringIn)
    if StringIn(i)=='.'
      SoundVec  = [SoundVec; sin(2*pi*SoundFreq*t); zeros(size(t))];
    end
    if StringIn(i)=='-'
      SoundVec  = [SoundVec; sin(2*pi*SoundFreq*t);sin(2*pi*SoundFreq*t); zeros(size(t))];
    end
    if StringIn(i)==' '
      SoundVec  = [SoundVec; zeros(size(t)); zeros(size(t))];
    end    
end

sound(SoundVec,fs)
        





end
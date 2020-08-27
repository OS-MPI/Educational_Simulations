function [M] = Langevin(B,C,varargin)

% This function is a basic magnetization function with external B field
% (H*mu0) in Tesla and the C(Concentration) as inputs. Optionally, you can
% specify the particle volume and saturation magnetization

Kb = 1.38e-23; %J/k

Ms = 446e3; %A/m
T = 300; %K
Vol = 4/3*pi*(12.5e-9)^3;%m^3 (12.5nm core radius)
B(B==0)=eps;



NumVarArgsIn = size(varargin,2);
if mod(NumVarArgsIn,2)>0
    error('Number of arguments inputted must be even')  
end

varargin = reshape(varargin,2,NumVarArgsIn/2);

for i = 1:NumVarArgsIn/2
    if strcmpi(varargin(1,i),'Ms')
        Ms = cell2mat(varargin(2,i));
    end
    
    if strcmpi(varargin(1,i),'Vol')
        Vol = cell2mat(varargin(2,i));
    end
    
    if strcmpi(varargin(1,i),'MoveSample') %Selecting MoveSample as 1(Yes) or 0(No) moves the sample along the axis to show how the signal changes with biasing
        MoveSample = cell2mat(varargin(2,i));
        
    elseif exist('MoveSampleSelected','var')
    else
        

    end
    
end




    M = C.*(coth(Vol*Ms*B./(Kb*T)) - (Kb*T)./(Vol*Ms*B));

    
end
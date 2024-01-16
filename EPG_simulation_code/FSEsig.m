function S = FSEsig(T2,B1,Gbg,opt)
%   EPG computation of the multiecho spin echo signal with shaped pulses
%   
%   Usage: S = FSEsig(T2,B1,Gbg,opt)
%   Author: R. Marc Lebel
%   Date: 07/2011
%   
%   Input:
%   T2: T2 relaxation time (s)
%   B1: scale factor for transmit (arbitrary units, near unity)
%   Gbg: Background field gradient (G/cm)
%   opt: options structure, defined in StimFit_optset and StimFit/StimFitImg
%       Must contain:
%           .mode: defines if selective or non-selective refocusing
%           .RFe : excitation pulse info (itself a structure)
%           .RFr : refocusing pulse info (itself a structure)
%           .etl : echo train length
%           .esp : echo spacing (s)
%           .Np  : number of points
%   
%   Output:
%   S: echo amplitudes (au)

%   Check inputs
if nargin ~= 4
    error('function requires 4 inputs');
end

%   Variable to replicate flip angles
replicate = ones(opt.etl,1);

%   Compute multi-echo signal
switch lower(opt.mode(1))
    
    %   Non-selective refocusing or 3D
    case 'n'
        
        %   Define vector of refocusing angles
        FA = B1*pi/180*opt.RFr.angle;
        FA = replicate*FA;
        
        %   Call compiled EPG code to compute echo train amplitudes
        S = epgMEX(opt.T1,T2,opt.esp,FA);
        
    %   2D mode with slice-selective refocus
    case 's'
        
        %   Compute magnetization following excitation
        Me = sin(B1*opt.RFe.alpha);
        
        %   Create array of flip angles at each echo train [etl x length(alpha)]
        %   This is assuming we're applying the same RF pulse throughout the train
        %   Modify if variable flip angles or different pulse shapes
        FA = replicate*(B1*opt.RFr.alpha);
        
        %   Call compiled EPG code to compute echo train amplitudes
        Mr = epgMEX(opt.T1,T2,opt.esp,FA);
        
        %   Combine excitation and refocusing
        %   NB: the relative decay during refocusing is independent of amount of
        %   excited magnetization, thus the excitation profile simply scales the
        %   echo amplitudes
        Mnet = (replicate*Me) .* Mr;
        
        %   Integrate signal across the slice profile and normalize
        %   (pretty arbitrary normalization...)
        S = sum(Mnet,2);
        S = S./opt.Nz;
        
end


%   Plot, if desired
if opt.debug
    
    %   Plotting axes
    te = opt.esp:opt.esp:opt.esp*opt.etl;
    
    %   Determine imaging mode 
    switch lower(opt.mode)
        case 'n'
            
            %   Create figure
            if ~any(findobj == 42575)
                figure(42575);
                scrsz = get(0,'ScreenSize');
                
                set(42575,'Name','EPG Simulated Magnetization','NumberTitle','off',...
                    'Position',[scrsz(3)/10 scrsz(4)/2 scrsz(3)/3 scrsz(4)/4],...
                    'Resize','off','Toolbar','none','Color','w');
            end
            
        case 's'
            
            %   Create figure
            if ~any(findobj == 42575)
                figure(42575);
                scrsz = get(0,'ScreenSize');
                set(42575,'Name','EPG Simulated Magnetization','NumberTitle','off',...
                    'Position',[scrsz(3)/10 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/2.25],...
                    'Resize','off','Toolbar','none','Color','w');
            end
            
            %   Plot excitation profile
            z = opt.Dz(1):(opt.Dz(2)-opt.Dz(1))/(opt.Nz-1):opt.Dz(2);
            subplot(2,4,1)
            plot(z,Me);
            set(gca,'XLim',opt.Dz,'YLim',[0 1],'LineWidth',1);grid on;
            xlabel('Position (cm)');ylabel('Mt (au)');title('Excitation Profile');
            
            %   Plot refocusing profiles (assuming uniform input magnetization)
            subplot(2,4,5)
            plot(z,Mr.');
            set(gca,'XLim',opt.Dz,'YLim',[0 1],'LineWidth',1);grid on;
            xlabel('Position (cm)');ylabel('Mt (au)');title('Refocusing Profiles');
            
            %   3D plot of magnetization decay
            %   Project the net signal onto an axis
            subplot(2,4,[2 3 6 7])
            [ETL,ALPHA] = ndgrid(te,z);
            surf(ETL,ALPHA,Mnet);grid on;box on;axis square;
            [ETL,ALPHA] = ndgrid(te,[0.97*opt.Dz(2) opt.Dz(2)]);
            Stmp = S .* mean(Mnet(:,round(opt.Nz/2))./S(:));
            hold on;surf(ETL,ALPHA,repmat(Stmp,[1 2]),'LineStyle','none');grid on;box on;hold off;
            set(gca,'YLim',opt.Dz,'XLim',[opt.esp opt.esp*opt.etl],'ZLim',[0 1],'LineWidth',1);
            view(50,25);
            ylabel('Position (cm)');xlabel('Echo time (s)');zlabel('Mt (au)');
            title('Combined Decay');
    end
    drawnow;
end

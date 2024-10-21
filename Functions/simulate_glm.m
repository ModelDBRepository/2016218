function [y stimcurr hcurr r] = simulate_glm(x,dt,k,h,dc,runs,SpikeTrains,BasisTab2,softRect,ind_plot)
% [y stimcurr hcurr r] = simulate_glm(x,dt,k,h,dc,runs)
%
%  --- Code Modified from Weber & Pillow (2017) ---:
%  Weber AI, Pillow JW. Capturing the Dynamical Repertoire of Single Neurons 
%  with Generalized Linear Models. Neural Comput. 2017 (12):3260-3289.

%  This code fits a Poisson GLM to given data, using basis vectors to
%  characterize the stimulus and post-spike filters.
%
%  The inputs are:
%   x: stimulus
%   dt: time step of x and y in ms
%   k: stimulus filter
%   h: post-spike filter
%   dc: dc offset
%   runs: number of trials to simulate
%   softRect: 0 uses exponential nonlinearity; 1 uses soft-rectifying nonlinearity
%   ind_plot: 0 for plotting, 1 for no plotting
%
%  The outputs are:
%   y: spike train (0s and 1s)
%   stimcurr: output of stimulus filter (without DC current added)
%   hcurr: output of post-spike filter
%   r: firing rate (stimcurr + hcurr + dc passed through nonlinearity)

%% generate data with fitted GLM
spikebinary = SpikeTrains;
nTimePts = length(x);
refreshRate = 1000/dt; % stimulus in ms, sampled at dt
if softRect
    NL = @logexp1;
else
    NL = @exp;
end

g = zeros(nTimePts+length(h),runs);     % filtered stimulus + dc
y = zeros(nTimePts,runs);               % initialize response vector (pad with zeros in order to convolve with post-spike filter)
r = zeros(nTimePts+length(h)-1,runs);   % firing rate (output of nonlinearity)
hcurr = zeros(size(g));     % post-spike current

stimcurr = sameconv(x,k);   % convolving new stimulus x, with fitted k filter
noise = 0 + 0.1.*randn(length(stimcurr),1);  % inject noise into filtered stimulus
Iinj = stimcurr + dc ; %+ noise;       % injected current includes DC drive

for runNum = 1:runs
    
    g(:,runNum) = [Iinj; zeros(length(h),1)]; 
    
    %%% loop to get responses, incorporate post-spike filter
    for t = 1:nTimePts
        r(t,runNum) = feval(NL,g(t,runNum));  % firing rate (output of nonlinearity)
        prob_0(t,runNum) = exp(-r(t,runNum)/refreshRate); % P(0 spikes) --- probability of 0 spikes in this time bin
        bern_prob(t,runNum) = 1-prob_0(t,runNum); % 1 - P(0 spikes) --- bernoulli 
        
%       % noise generation (random walk -- gaussian centered at 0)
%         if t == 1
%             rand_num(t,runNum) = 0 + 0.005*randn;
%         elseif t>1
%             rand_num(t,runNum) = rand_num(t-1,runNum) + (0 + 0.005*randn);
%         end
        
        % spike generation: 
        if 0.001+rand<bern_prob(t,runNum)     % if probability of spiking is greater than randomly drawn number (ADDED 0.001 SHIFT TO GET RID OF NOISE!!)
           y(t,runNum) = 1;
           g(t:t+length(h)-1,runNum) = g(t:t+length(h)-1,runNum) + h;  % add post-spike filter
           hcurr(t:t+length(h)-1,runNum) = hcurr(t:t+length(h)-1,runNum) + h;
           %rand_num(t,runNum) = 0; % reset random number walk back to 0
        end
    end
end

hcurr = hcurr(1:nTimePts,:);  % trim zero padding for post-spike current
r = r(1:nTimePts,:);  % trim zero padding for firing rate

%% plot the results 
% % time
minT = 1/dt;
maxT = length(x);
tIdx = minT:maxT;
t = (tIdx-minT)*dt;
% % 
% figure('Position', [10 10 1200 650])
% %%% stimulus
% yo(1) = subplot(5,1,1); hold on;
% plot(t,x(tIdx),'color','k','linewidth',2)
% xlim([min(t) max(t)])
% ylim([min(x(tIdx))-.05*abs(min(x(tIdx))) max(x(tIdx))+.05*abs(max(x(tIdx)))])
% box off
% title('stimulus')
% 
% %%% filter outputs
% yo(2) = subplot(5,1,2); hold on;
% plot(t,Iinj(tIdx),'r','linewidth',1.5);   % plotting the stimulus current (+dc)
% plot(t,hcurr(tIdx,1),'b','linewidth',1.5) % plotting the post-spike current
% plot(t,g(tIdx,1),'g','linewidth',1.5)
% xlim([min(t) max(t)])
% ylim([min([hcurr(tIdx,1); Iinj(tIdx)])*1.1 max([hcurr(tIdx,1); Iinj(tIdx)])*1.1])
% box off
% title('filter outputs')
% 
% %%% firing rate (lambda)
% yo(3) = subplot(5,1,3); hold on;
% semilogy(t,feval(NL,hcurr(tIdx,1)+Iinj(tIdx)),'color',[.5 .5 .5],'linewidth',1.5)
% xlim([min(t) max(t)])
% box off
% title('lambda (conditional intensity/IFR)')
% 
% %%% probability of spiking (between 0 and 1)
% yo(6)=subplot(5,1,4);
% plot(t,prob_0(tIdx),'color',[0.5,0.5,0.5]);
% title('P(spike|lambda) = 1 - exp(-lambda/refresh rate)')
% ylim([0,1])
% xlim([min(t) max(t)])
% 
% %%% GLM spikes
% yo(4) = subplot(5,1,5); hold on;
% spikeHeight = .7;
% 
% for i = 1:size(y,2) % for each run of glm simulation
%     spt = find(y(tIdx,i));
%     for spikeNum = 1:length(spt)
%         plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i-.5 i-.5+spikeHeight],'color',[.5 .5 .5],'linewidth',1.25)
%     end
% end
% 
% for f = 1:size(spikebinary,2)   % for each real spike
%     spt = find(spikebinary(tIdx,f));
%     for spikeNum = 1:length(spt)
%         plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i+f+.5 i+f+.5+spikeHeight],'color',[0 0 0],'linewidth',1.25)
%     end
% end
% 
% xlim([0 max(t)-min(t)])
% ylim([0 runs+spikeHeight+2])
% xlabel('time (ms)')
% title('spikes -- if P>rand, a spike occurs')
% linkaxes(yo,'x')
   

%% Second Figure:
%%% plot filters and compare real to model spikes
if ind_plot == 0
    axisLabelFontSize = 12;
    axisTickLabelFontSize = 12;
    axisWidth = 1;
    
    fig=figure; fig.Position=[10 10 1200 650];
    subplot(3,4,8); hold on;
    h1 = gca;
    plot([0 length(k)],[0 0],'k--','linewidth',1.5);
    plot(k,'b','linewidth',2);
    set(gca,'xtick',0:length(k)/4:length(k),'xticklabel',round(-length(k)*dt:length(k)/4*dt:0))
    set(gca,'tickdir','out','linewidth',axisWidth,'fontsize',axisTickLabelFontSize)
    text(length(k)/15,max(k)-.05*(max(k)-min(k)),['\mu = ' num2str(round(dc*10)/10)],'fontsize',axisLabelFontSize);
    xlim([-5 length(k)])
    ylim([min(k)-.05*(max(k)-min(k)) max(k)+.05*(max(k)-min(k))])
    box off
    h1p = get(h1,'position');
    %set(h1,'position',[h1p(1) h1p(2)*1.05 h1p(3)*.9 h1p(4)*.95])
    
    subplot(3,4,12); hold on;
    plot([0 length(h)],[0 0],'k--','linewidth',1.5);
    plot(h(2:end),'r','linewidth',2);
    h2 = gca;
    set(gca,'tickdir','out','xtick',0:length(h)/4:length(h),'xticklabel',round(0:(length(h)/4*dt):length(h)*dt))
    set(gca,'linewidth',axisWidth,'fontsize',axisTickLabelFontSize)
    xlabel('time (ms)','fontsize',axisLabelFontSize)
    xlim([-5 length(h)])
    ylim([min(h)-.05*(max(h)-min(h)) max(h)+.05*(max(h)-min(h))])
    box off
    h2p = get(h2,'position');
    %set(h2,'position',[h2p(1) h2p(2)*.95 h2p(3)*.9 h2p(4)*.95])
    
    fig(1)=subplot(3,4,1:3); 
    plot(t,x(tIdx),'color','k','linewidth',2)
    %xlim([55400 56600])
    xlim([min(t) max(t)])
    ylim([min(x(tIdx))-.05*abs(min(x(tIdx))) max(x(tIdx))+.05*abs(max(x(tIdx)))])
    box off
    title('stimulus')
    
    fig(2)=subplot(3,4,[5:7,9:11]); hold on;
    spikeHeight = .7;
    
    for i = 1:size(y,2) % for each run of glm simulation
        spt = find(y(tIdx,i));
        for spikeNum = 1:length(spt)
            plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i-.5 i-.5+spikeHeight],'color',[.5 .5 .5],'linewidth',1.25)
        end
    end
    
    for f = 1:size(spikebinary,2)   % for each real spike
        spt = find(spikebinary(tIdx,f));
        for spikeNum = 1:length(spt)
            plot([spt(spikeNum)*dt spt(spikeNum)*dt],[i+f-.5 i+f-.5+spikeHeight],'color',[0 0 0],'linewidth',1.25)
        end
    end
    
    ylim([0 runs+spikeHeight+1])
    %xlim([55400 56600])
    xlim([0 max(t)])
    xlabel('time (ms)')
    title('spikes (black=real, grey=GLM)')
    linkaxes(fig,'x')
    
    uit = uitable('Data', table2cell(BasisTab2),'ColumnName',BasisTab2.Properties.VariableNames,...
        'Units', 'Normalized', 'Position',[0.72,0.75,0.2,0.15]);
    set(uit,'ColumnWidth',{50},'FontSize',13)
end



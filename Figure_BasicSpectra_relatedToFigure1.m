% Required inputs :
% Wake, SWSEpoch, REMEpoch, Epoch: intervalSets extracted from the sleep scoring (see SleepScoring_Accelero_OBGAmma.m)
% Spectrum of the region of interest (.mat)

 
clear all
Cols = {'k','r'};


SpeedBins = [2:2:30];

for grp = 1:2 %control mice and social defeat stress (SDS) mice 
    for mm = 1:length(enregistrements{grp}) % enregistrements : file path for each animal
        cd(enregistrements{grp}{mm})
        mm
        % load SeepScoring
        load('SleepScoring_Accelero.mat','REMEpoch','Wake','SWSEpoch')
        % Restrict to the same duration for all mice
        Epoch = intervalSet(0*3600*1e4,7.5*3600*1e4);
        Wake = and(Wake,Epoch);
        REMEpoch = and(REMEpoch,Epoch);
        SWSEpoch = and(SWSEpoch,Epoch);
        
        load('behavResources.mat', 'Vtsd') %speed variable (tsd format) extracted from video tracking
        SlowWake = and(Wake,thresholdIntervals(Vtsd,4,'Direction','Below'));
        FastWake = and(Wake,thresholdIntervals(Vtsd,4,'Direction','Above'));
        
        clear Spectro
        load('Bulb_deep_Low_Spectrum.mat') % Olfactory bulb spectrogram extracted from the sleepscoring
        Sptsd = tsd(Spectro{2}*1e4,(Spectro{1}));
        MeanSpec.OB.SWS(grp,mm,:) = nanmean(Data(Restrict(Sptsd,SWSEpoch)));
        MeanSpec.OB.REM(grp,mm,:) = nanmean(Data(Restrict(Sptsd,REMEpoch)));
        MeanSpec.OB.Wake(grp,mm,:) = nanmean(Data(Restrict(Sptsd,Wake)));
        
        % Speed Matched
        clear Spectemp
        for sp = 1:length(SpeedBins)-1
            LittleWake = and(and(Wake,thresholdIntervals(Vtsd,SpeedBins(sp+1),'Direction','Below')),thresholdIntervals(Vtsd,SpeedBins(sp),'Direction','Above'));
            Spectemp(sp,:)= nanmean(Data(Restrict(Sptsd,LittleWake)));
        end
        MeanSpec.OB.FastWake(grp,mm,:) = nanmean(Spectemp,1);
        
        f{1} = Spectro{3};
        
        clear Spectro
        if exist('dHPC_rip_Low_Spectrum.mat') % HPC spectrogram extracted from the sleepscoring
            disp('HPC rip')
            load('dHPC_rip_Low_Spectrum.mat')
            Sptsd = tsd(Spectro{2}*1e4,(Spectro{1}));
            MeanSpec.dHPC.SWS(grp,mm,:) = nanmean(Data(Restrict(Sptsd,SWSEpoch)));
            MeanSpec.dHPC.REM(grp,mm,:) = nanmean(Data(Restrict(Sptsd,REMEpoch)));
            MeanSpec.dHPC.Wake(grp,mm,:) = nanmean(Data(Restrict(Sptsd,Wake)));
            clear Spectemp
            for sp = 1:length(SpeedBins)-1
                LittleWake = and(and(Wake,thresholdIntervals(Vtsd,SpeedBins(sp+1),'Direction','Below')),thresholdIntervals(Vtsd,SpeedBins(sp),'Direction','Above'));
                Spectemp(sp,:)= nanmean(Data(Restrict(Sptsd,LittleWake)));
            end
            MeanSpec.dHPC.FastWake(grp,mm,:) = nanmean(Spectemp,1);
            
            f{2} = Spectro{3};
        else
            MeanSpec.dHPC.SWS(grp,mm,:) = zeros(1,261);
            MeanSpec.dHPC.REM(grp,mm,:) = zeros(1,261);
            MeanSpec.dHPC.Wake(grp,mm,:) = zeros(1,261);
            MeanSpec.dHPC.FastWake(grp,mm,:) = zeros(1,261);
            
        end
        
        clear Spectro
        try,load('PFCx_sup_Low_Spectrum.mat') %PFC spectrogram 
        catch
            load('PFCx_deep_Low_Spectrum.mat')
        end
        f{3} = Spectro{3};
        Sptsd = tsd(Spectro{2}*1e4,(Spectro{1}));
        MeanSpec.PFC.SWS(grp,mm,:) = nanmean(Data(Restrict(Sptsd,SWSEpoch)));
        MeanSpec.PFC.REM(grp,mm,:) = nanmean(Data(Restrict(Sptsd,REMEpoch)));
        MeanSpec.PFC.Wake(grp,mm,:) = nanmean(Data(Restrict(Sptsd,Wake)));
        clear Spectemp
        for sp = 1:length(SpeedBins)-1
            LittleWake = and(and(Wake,thresholdIntervals(Vtsd,SpeedBins(sp+1),'Direction','Below')),thresholdIntervals(Vtsd,SpeedBins(sp),'Direction','Above'));
            Spectemp(sp,:)= nanmean(Data(Restrict(Sptsd,LittleWake)));
        end
        MeanSpec.PFC.FastWake(grp,mm,:) = nanmean(Spectemp,1);    end
end

close all

Regions = fieldnames(MeanSpec); %OB / HPC / PFC
Epochs = fieldnames(MeanSpec.PFC); %SWS / REM / Wake / Fast Wake

for ep = 1:length(Epochs)
    fig = figure;
    for reg = 1:length(Regions)
        for IndivMice = 0:1
            subplot(2,3,reg+3*IndivMice)
            % Ctrl
            dat = squeeze(MeanSpec.(Regions{reg}).(Epochs{ep})(1,:,:));
            dat(find(sum(dat')==0),:) = [];
            dat =  dat./repmat(mean(dat(:,50:end)'),size(dat,2),1)';
            if IndivMice ==1
                plot(f{reg},dat,'k')
            else
                shadedErrorBar(f{reg},nanmean(log(dat)),(stdError(log(dat))),'k')
            end
            hold on
            MnDel{reg}{1} = max(dat(:,1:find(f{reg}>3,1,'first'))');
            MnThet{reg}{1} = max(dat(:,find(f{reg}>5,1,'first'):find(f{reg}>10,1,'first'))');
            
            % SDS
            dat = squeeze(MeanSpec.(Regions{reg}).(Epochs{ep})(2,:,:));
            dat(find(sum((dat==0)')),:) = [];
            dat =  dat./repmat(mean(dat(:,50:end)'),size(dat,2),1)';
            if reg ==3
                % This mouse has a huge 1/f{reg} so excluded as outlier
                dat(end,:) = [];
            end
            if IndivMice ==1
                plot(f{reg},dat,'r')
            else
                shadedErrorBar(f{reg}, nanmean(log(dat)),(stdError(log(dat))),'r')
            end
            MnDel{reg}{2} = max(dat(:,1:find(f{reg}>3,1,'first'))');
            MnThet{reg}{2} = max(dat(:,find(f{reg}>5,1,'first'):find(f{reg}>10,1,'first'))');
            
            title(Regions{reg})
            axis square
            xlabel('Frequency')
            ylabel('Power')
            makepretty
        end
    end
    %         saveas(fig.Number,['Spectra',Epochs{ep},'.fig'])
    
    fig = figure;
    for reg = 1:length(Regions)
        subplot(2,3,reg)
        MakeSpreadAndBoxPlot2_SB(MnDel{reg},{},[1:2],{'Ctrl','SDS'},'paired',0)
        title(Regions{reg})
        ylabel('DeltaPower')
        makepretty
        subplot(2,3,reg+3)
        MakeSpreadAndBoxPlot2_SB(MnThet{reg},{},[1:2],{'Ctrl','SDS'},'paired',0)
        title(Regions{reg})
        ylabel('ThetaPower')
        makepretty
    end
    %     saveas(fig.Number,['SpecQuant',Epochs{ep},'.fig'])
end
%% read kilosort to singleunits
% cgs - 0 = noise / 1 = mua / 2 = good / 3 = unsorted
clear;
addpath(genpath('D:\Kilosort2-user\my_kilosort'));
Datapath={
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\Intan\';...
    };
Laser_artifact={'n';...
    'n';...
    'n';...
    'n';...
    };

for k=1:length(Datapath)
    laser_artifact=Laser_artifact{k};
    for probe=1:2
        for shank=1:2
            istr = strfind(Datapath{k},'\');
            datapath=[Datapath{k},'Kilosort_data',num2str(probe),num2str(shank),'\'];
            % 在session日期路径下保存single units
            [suapath_root, ~, ~] = fileparts(Datapath{k});
            [suapath_root, ~, ~] = fileparts(suapath_root);
            suapath = [suapath_root, '\Kilosort_SingleUnits', num2str(probe), '\'];  % shank体现在后文文件名中
            if ~exist(suapath, 'dir')
                mkdir(suapath);
            end
            rawdatafile=dir([datapath,'/*.bin']);
            trialtimefile=dir([datapath,'Trial_time_all*.mat']);%_new
            load([datapath,trialtimefile.name]);
            spike_sample_time = readNPY([datapath,'spike_times.npy']);
            spike_clusters = readNPY([datapath,'spike_clusters.npy']);
            spike_amplitudes = readNPY([datapath,'amplitudes.npy']);
            [cidstmp, cgs] = readClusterGroupsCSV([datapath,'cluster_group.tsv']);  %mua 1; good 2; noise 3
            [cids, channel, firing_rate]= readClusterInfoCSV([datapath,'cluster_info.tsv']);  % already +1 to be 1-384, 0-383 in kilosort phy, cids is from min to max order
            chMap = readNPY(fullfile(datapath, 'channel_map.npy'))+1;  % Order in which data was streamed to disk; must be 1-indexed for Matlab
            chPos = readNPY(fullfile(datapath, 'channel_positions.npy')); 

            %keep only good units
            cids_good=cidstmp(find(cgs==2));
            if ~isempty(cids_good)
                channel_good=channel(find(ismember(cids,cids_good)));%already 1-384 range
                firing_rate_good=firing_rate(find(ismember(cids,cids_good)));
                idx=find(ismember(spike_clusters,cids_good));
                spike_clusters_good=spike_clusters(idx);
                spike_sample_time_good=spike_sample_time(idx);
                spike_amplitudes_good=spike_amplitudes(idx);
                
                gwfparams.dataDir = datapath;    % KiloSort/Phy output folder
                gwfparams.fileName = rawdatafile.name;         % .dat file containing the raw 
                gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)

                if max(chMap)<=32
                % CB64=2x32
                    fs=20000;
                    gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
                    gwfparams.wfWin = [-10 29];              % Number of samples before and after spiketime to include in waveform %[-10 29] for CB64  [-15 44] for NP 30K fs
                    gwfparams.nWf = 10000;                    % Number of waveforms per unit to pull out
                else
                    % NP384
                    fs=30000;
                    gwfparams.nCh = 384;                     % Number of channels that were streamed to disk in .dat file
                    gwfparams.wfWin = [-15 44];              % Number of samples before and after spiketime to include in waveform %[-10 29] for CB64  [-15 44] for NP 30K fs
                    gwfparams.nWf =1000;                    % Number of waveforms per unit to pull out
                end

                gwfparams.spikeTimes = spike_sample_time_good; % Vector of cluster spike times (in samples) same length as .spikeClusters
                gwfparams.spikeClusters = spike_clusters_good; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
                
                wf = getWaveForms(gwfparams);
                
                figure;
                for j=1:length(wf.unitIDs)
                    chj=channel_good(find(cids_good==wf.unitIDs(j)));%already 1-384 range
                    unit.ID=wf.unitIDs(j);
                    unit.firing_rate=firing_rate_good(find(cids_good==wf.unitIDs(j)));
                    chjtmp=find(chMap==chj);
                    unit.waveforms=squeeze(wf.waveForms(j,:,chjtmp,:));
                    Distance=sqrt(sum((chPos-repmat(chPos(chjtmp,:),size(chPos,1),1)).^2,2));% update by Guang 20250329
                    [ss,ii]=sort(Distance);
                    unit.waveforms_allch=squeeze(wf.waveForms(j,:,ii(1:32),:));
                    unit.wvsch=chMap(ii(1:32));
                    unit.wvspos=chPos(ii(1:32),:);

                    idx_tmp=find(spike_clusters_good==wf.unitIDs(j));
                    idx_tmp=idx_tmp(1:end-1);
                    unit.spike_times=Trial_time_all(1,spike_sample_time_good(idx_tmp))';
                    unit.trials=Trial_time_all(2,spike_sample_time_good(idx_tmp))';
                    unit.spike_times_wv_keeps=Trial_time_all(1,wf.spikeTimeKeeps(j,find(~isnan(wf.spikeTimeKeeps(j,:)))));
                    unit.trials_wv_keeps=Trial_time_all(2,wf.spikeTimeKeeps(j,find(~isnan(wf.spikeTimeKeeps(j,:)))));
                    unit.amplitudes=spike_amplitudes_good(idx_tmp);

                    if length(unit.trials)~=length(unit.amplitudes)
                       'error: spike amp mismatch';
                    end
                    if shank==1
                       unit.channel=repmat(chj,length(unit.spike_times),1);%DBC 32*2 shank1  or NP 1.0 or 2.0 %already 1-384 range
                    else
                       unit.channel=repmat(chj+32,length(unit.spike_times),1);%DBC 32*2 shank2
                    end
                    unit.stable_trials=min(unit.trials):max(unit.trials);
                    unit = func_classify_unit_intan(unit,fs,0);
                    unit.manual_quality_score='1';
                    unit.artifact_annotation=laser_artifact;
                    if exist('oepxi_motor12')
                        unit.oepxi_motor12=oepxi_motor12;
                    else
                        unit.oepxi_motor12=[];
                    end

                    disp(['saving ',suapath, ' unit ',num2str(j)]);
                    save([suapath, 'SingleUnit',num2str(shank),'_',num2str(j),'.mat'],'unit');

                    if unit.cell_type==1
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'r');hold on;
                        subplot(2,4,2);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'r-'); hold on;
                    elseif unit.cell_type==2
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'b');hold on;
                        subplot(2,4,3);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'b-'); hold on;
                    else
                        subplot(2,4,1);hold on;
                        scatter(rand(1)*10,unit.wvspos(1,2)-max(chPos(:,2)),'k');hold on;
                        subplot(2,4,4);hold on;
                        plot(mean(unit.waveforms, 'omitnan'), 'k-'); hold on;
                    end
                    subplot(2,4,5:8);hold on;
                    scatter(unique(unit.trials),j*ones(1,length(unique(unit.trials))),'g+');hold on;
                    xlabel('trial');
                    ylabel('unit');

                end
                saveas(gcf, [suapath, 'celltype_depth_map.png'],'png');
                close;
                clearvars -except Datapath Laser_artifact laser_artifact k position probe shank
            end
        end
    end
end

%% manual correct Spike Timing kilosort
clear;
close all

session_dir = {
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\';...
     };

for probe=1:2
    for i_session = 1:size(session_dir,1)
        single_unit_dir = [session_dir{i_session},'Kilosort_SingleUnits',num2str(probe),'\'];
        unit_files = dir([single_unit_dir,'*.mat']);
        for i_unit = 1:size(unit_files,1)
            load([single_unit_dir, unit_files(i_unit).name]); 
            unit_tmp = unit;
            isi = diff(unit_tmp.spike_times);
            i_1st_spk = diff(unit_tmp.trials)>0;
           
            i_spk_discard = find(isi<.0005  & i_1st_spk==0);

            unit_tmp.spike_times(i_spk_discard,:) = [];
            unit_tmp.amplitudes(i_spk_discard,:) = [];
            unit_tmp.trials(i_spk_discard,:) = [];
            unit_tmp.channel(i_spk_discard,:) = [];
            unit_tmp.stable_trials   =  unit.stable_trials;
            unit_tmp.cell_type = unit.cell_type;
            unit_tmp.ID=unit.ID;
            unit_tmp.firing_rate=unit.firing_rate;
            unit_tmp.manual_quality_score=unit.manual_quality_score;
            unit_tmp.artifact_annotation=unit.artifact_annotation;

            disp([num2str(size(i_spk_discard,1)),'/',num2str(size(isi,1)+1),' spikes discarded'])
            unit = unit_tmp;
            ISI = diff(unit.spike_times);
            ISI = ISI(find(ISI<.5));
            ISI = [ISI; -ISI];
            unit.false_alarm_est = sum(abs(ISI)<.002)/length(ISI);
       
            unit.waveforms=unit.waveforms(find(~isnan(unit.waveforms(:,1))),:);
            wave_amp_tmp = range(unit.waveforms,2);
            if size(wave_amp_tmp,1)>100
                mean_wave_amp = conv(wave_amp_tmp,ones(1,100)/100,'same');
                mean_wave_amp(1:50) = mean_wave_amp(51);
                mean_wave_amp(end-49:end) = mean_wave_amp(end-50);
                wave_amp_tmp = wave_amp_tmp-mean_wave_amp;   % only look at the residues
            end
       
            mu_est = mean(wave_amp_tmp);
            sigma_est = std(wave_amp_tmp);
       
            X_min = (mu_est-sigma_est*5);
            X_max = (mu_est+sigma_est*5);
            X = X_min:(X_max-X_min)/100:X_max;
            Y_fit = normpdf(X,mu_est,sigma_est);
            Y = histcounts(wave_amp_tmp, X);  
            Y(end+1) = sum(wave_amp_tmp == X(end)); % histcounts 返回长度为 length(edges)-1，需要手动加上最后一个 bin
            Y = Y/sum(Y)/((X_max-X_min)/100);
            r = corr(Y,Y_fit');
       
            unit.miss_est = r^2;
            save([single_unit_dir, unit_files(i_unit).name],'unit'); 

        end
    end
clearvars -except session_dir position probe
end

%% manual check for duplicate clusters kilosort
clear;
close all

session_dir = {
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\';...
     };

probetype=1; %1 CB64; 2 NP1.0 or NP2.0，这是CB64记录的

for probe=1:2
    for i_session = 1:size(session_dir,1)
        single_unit_dir = [session_dir{i_session},'Kilosort_SingleUnits',num2str(probe),'\'];%num2str(position),
        if ~exist([single_unit_dir,'CheckDuplicate'])
            mkdir([single_unit_dir,'CheckDuplicate'])
        end

        unit_files = dir([single_unit_dir,'SingleUnit*.mat']);
        for i_unit = 1:size(unit_files,1)
            load([single_unit_dir, unit_files(i_unit).name]);
            unit_tmp = unit;
            unit_channel = median(unit_tmp.channel);
            unit1pos=unit.wvspos(1,:);

            for i_unit_pair = i_unit+1:size(unit_files,1)
                load([single_unit_dir, unit_files(i_unit_pair).name]);
                unit_pair_tmp = unit;
                unit_pair_channel = median(unit_pair_tmp.channel);
                unit2pos=unit.wvspos(1,:);

                if strcmp(unit_files(i_unit).name(11),unit_files(i_unit_pair).name(11))
                    if sqrt(sum((unit1pos-unit2pos).^2))<30 %um %abs(unit_channel-unit_pair_channel)<=1
                        if length(unique(unit_tmp.trials))>10&length(unique(unit_pair_tmp.trials))>10
                            close all
                            func_check_for_duplicate_unit_kilosort(unit_tmp, unit_files(i_unit).name(1:end-4), unit_pair_tmp, unit_files(i_unit_pair).name(1:end-4), single_unit_dir, probetype);
                        end
                    end
                end  
            end
        end
    end
    clearvars -except session_dir position probe probetype
end

%% mannually put in the duplicate unit pair (merge single unit for each subfolder)
clear;
single_unit_dir='C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\Kilosort_SingleUnits1\';
combine_pairs={
  % '1_11' '1_8';...
        };

if size(combine_pairs,1)>0    

    for i_pair = 1:size(combine_pairs,1)
        tf=0;
        pair1=combine_pairs{i_pair,1};
        pair2=combine_pairs{i_pair,2};
        if exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tf=1;
        elseif ~exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tmp=[];
            for k=1:size(combine_pairs,1)
                tmp(k)=strcmp(combine_pairs{k,2},pair1);
            end
            tmp=find(tmp);
            pair1=combine_pairs{tmp(1),1};
            tf=1;
        elseif exist([single_unit_dir, 'SingleUnit',pair1,'.mat'],'file')&&~exist([single_unit_dir, 'SingleUnit',pair2,'.mat'],'file')
            tmp=[];
            for k=1:size(combine_pairs,1)
                tmp(k)=strcmp(combine_pairs{k,2},pair2);
            end
            tmp=find(tmp);
           pair1=combine_pairs{tmp(1),1};
           tf=1;
        end
        if tf
            disp(['combining SingleUnit ',pair1,'  ',pair2]);
            load([single_unit_dir, 'SingleUnit',pair1,'.mat']);
            unit1 = unit;
            load([single_unit_dir, 'SingleUnit',pair2,'.mat']);
            unit2 = unit;
            clear unit;
            [unit] = func_combine_duplicate_unit_kilosort(unit1,unit2);
            isi = diff(unit.spike_times);
            i_1st_spk = diff(unit.trials)>0;
            i_spk_discard = find(isi<.0005  & i_1st_spk==0);
            
            unit.spike_times(i_spk_discard,:) = [];
            unit.amplitudes(i_spk_discard,:) = [];
            unit.trials(i_spk_discard,:) = [];
            unit.channel(i_spk_discard,:) = [];
            
            disp([num2str(size(i_spk_discard,1)),'/',num2str(size(isi,1)+1),' spikes discarded'])
            
            % score unit quality
            ISI = diff(unit.spike_times);
            ISI = ISI(ISI<.5);
            ISI = [ISI; -ISI];
            unit.false_alarm_est = sum(abs(ISI)<.002)/length(ISI);
            
            unit.waveforms=unit.waveforms(find(~isnan(unit.waveforms(:,1))),:);
            wave_amp_tmp = range(unit.waveforms,2);
            if size(wave_amp_tmp,1)>100
                mean_wave_amp = conv(wave_amp_tmp,ones(1,100)/100,'same');
                mean_wave_amp(1:50) = mean_wave_amp(51);
                mean_wave_amp(end-49:end) = mean_wave_amp(end-50);
                wave_amp_tmp = wave_amp_tmp-mean_wave_amp;   % only look at the residues
            end
            
            mu_est = mean(wave_amp_tmp);
            sigma_est = std(wave_amp_tmp);
            
            X_min = (mu_est-sigma_est*5);
            X_max = (mu_est+sigma_est*5);
            
            X = X_min:(X_max-X_min)/100:X_max;
            Y_fit = normpdf(X,mu_est,sigma_est);
            Y = histcounts(wave_amp_tmp,X);
            Y(end+1) = sum(wave_amp_tmp == X(end));
            Y = Y/sum(Y)/((X_max-X_min)/100);
            r = corr(Y,Y_fit');
            
            unit.miss_est = r^2;
            
            disp(['Updateing SingleUnit', pair1]);
            save([single_unit_dir, 'SingleUnit',pair1,'.mat'], 'unit');
    
            disp(['Deleting SingleUnit', pair2]);
            delete([single_unit_dir, 'SingleUnit',pair2,'.mat']);
            if exist([single_unit_dir, 'SingleUnit',pair2,'.png'],'file')
               delete([single_unit_dir, 'SingleUnit',pair2,'.png']);
            end
        end
    end
end

%% manual final plot and check badunits kilosort
clear;
close all

session_dir = {
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\';...
     };

for probe=1:2
    for i_session = 1:size(session_dir,1)
        single_unit_dir = [session_dir{i_session},'Kilosort_SingleUnits',num2str(probe),'\'];
        if exist(single_unit_dir)
            if ~exist([single_unit_dir,'BadUnits'])
                mkdir([single_unit_dir,'BadUnits'])
            end
            unit_files = dir([single_unit_dir,'*.mat']);
    
            for i_unit = 1:size(unit_files,1)
                load([single_unit_dir, unit_files(i_unit).name]); 
                figure;
                subplot(2,3,1);
                isi = diff(unit.spike_times);
                isi = isi(find(isi<.5));
                isi = [isi; -isi];
                edges = -.03:.00025:.03;
                n = histcounts(isi, edges);       % 长度为 length(edges)-1
                n(end+1) = sum(isi == edges(end)); 
                plot(edges, n, 'r')
                if max(n)~=0
                    axis([-.02 .02 0 max(n)]);
                end
                xlabel('s');
                ylabel('count');
                subplot(2,3,4);
                if unit.cell_type == 1
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'r')
                end
                if unit.cell_type == 2
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'b')
                end
                if unit.cell_type == 0
                    plot(mean(unit.waveforms, 1, 'omitnan'), 'k')
                end
                xlabel('sample count');
                ylabel('uV');
       
                spike_times_psth = {};
                n_trial = 0;
                trial0=min(unit.trials):max(unit.trials);
                trial_tmp=intersect(1:999,trial0)'; % only behavior trial

                if size(trial_tmp,1)>1
                    trial_tmp=trial_tmp';
                end
                if ~isempty(trial_tmp)
                    for i_trial = trial_tmp
                        n_trial = n_trial+1;
                        spike_times_psth{n_trial,1} = unit.spike_times(unit.trials==i_trial)';
                    end
                end
                if size(spike_times_psth,1)>10&max(unit.spike_times)>0
                    subplot(4,3,5);
                    [psth, t] = func_getPSTH(spike_times_psth,0,max(unit.spike_times));
                    bar(t,psth,'k');hold on;
                    xlim([0 min([10 max(unit.spike_times)])])
                    xlabel('s');
                    ylabel('spikes/s');
                    subplot(4,3,2);
                    trial_idx=find(ismember(unit.trials,trial_tmp));
                    plot(unit.spike_times(trial_idx),unit.trials(trial_idx),'.k');
                    xlim([0 min([10 max(unit.spike_times)])])
                    xlabel('s');
                    ylabel('trial');
                end
                subplot(4,3,[3 6]);
                trial_idx=find(unit.trials>=1000);
                plot(unit.spike_times(trial_idx),unit.trials(trial_idx),'.k');hold on;
                xlim([-10 10])
                xlabel('s');
                ylabel('trial');

                saveas(gcf,[single_unit_dir, unit_files(i_unit).name(1:end-4),'.png'],'png');
                close;
            end
        end
    end
end

%% remove bad bitcode trial

clear;
close all
session_dir = {
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\';...
     };

for probe=1:2
    for i_session = 1:size(session_dir,1)
    single_unit_dir = [session_dir{i_session},'Kilosort_SingleUnits',num2str(probe),'\'];
    unit_files = dir([single_unit_dir,'*.mat']);
        for i_unit = 1:size(unit_files,1)
           load([single_unit_dir, unit_files(i_unit).name]); 
           idx=unique(unit.trials);
           idx=find(unit.trials<=idx(5)|unit.trials>=idx(end-4));
           unit.spike_times(idx)=[];
           unit.trials(idx)=[];
           unit.amplitudes(idx)=[];
           unit.channel(idx)=[];
           unit.stable_trials=min(unit.trials):max(unit.trials);
           save([single_unit_dir, unit_files(i_unit).name],'unit');
        end
    end
end

%% compile obj first without unit
aom=max(obj.wavesurfer.aom_input_trace,[],2);
x=max(obj.wavesurfer.xGalvo_trace,[],2);
y=max(obj.wavesurfer.yGalvo_trace,[],2);
SCsti_trial=find(aom>0.2&x<0.2&y<0.2);
SNrsti_trial=find(aom>0.2&x>0.2&y<0.2);
ALMsti_trial=find(aom>0.2&x<0.2&y>0.2);
control_trial=find(aom<=0.2);
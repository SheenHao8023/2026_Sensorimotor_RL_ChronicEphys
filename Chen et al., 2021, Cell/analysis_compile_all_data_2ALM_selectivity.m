% establish, compute and save trial type selectivity variables
clear all
close all

addpath('.\helper_functions\')
filepath='C:\Users\XinHao\Desktop';

% list of animals to analyze
Animals_list = {
    % standard task
    'BAYLORGC4';...  % VGAT
    'BAYLORGC12';...  % PVReachR
    'BAYLORGC13';...  % PVReachR
    'BAYLORGC15';...  % PVReachR
    'BAYLORGC17';...  % PVReachR    
    'BAYLORGC18';...  % PVReachR
    'BAYLORGC19';...  % PVReachR
    'BAYLORGC21';...  % PVReachR
    'BAYLORGC22';...  % PVReachR
    'BAYLORGC25';...  % PVReachR
    'BAYLORGC26';...  % PVReachR
    'BAYLORGC50';...  % VGAT
    'BAYLORGC72';...  % VGAT
    'BAYLORGC75';...  % VGAT
    'BAYLORGC86';...  % PVReachR
    'BAYLORGC87';...  % PVReachR
    'BAYLORGC93';...  % PVReachR
    'BAYLORGC95';...  % PVReachR
    
    % reversed contingency task
    'BAYLORGC61';...  % VGAT
    'BAYLORGC69';...  % VGAT
    'BAYLORGC71';...  % VGAT
    'BAYLORGC74';...  % VGAT
    'BAYLORGC78';...  % VGAT
    
    % fully reversed task
    'BAYLORGC79';...  % PVReachR
    'BAYLORGC80';...  % PVReachR
    'BAYLORGC84';...  % PVReachR
    'BAYLORGC85';...  % PVReachR
    'BAYLORGC83';...  % PVReachR
    'BAYLORGC88';...  % PVReachR
    'BAYLORGC81';...  % PVReachR
    'BAYLORGC89';...  % PVReachR
    'BAYLORGC90';...  % PVReachR
    'BAYLORGC94';...  % PVReachR
    
    % reversed tactile stimulus task
    'BAYLORGC91';...  % PVReachR
    'BAYLORGC92';...  % PVReachR
    'BAYLORGC96';...  % PVReachR
    'BAYLORGC98';...  % PVReachR
    'BAYLORGC99';...  % PVReachR
    'BAYLORGC100';...  % PVReachR
    
    };
N_animals = length(Animals_list);  

for I_animal = 1:N_animals
    obj=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_PSTH_',Animals_list{I_animal},'.mat'],...
                    'N_trials_all','Sig_selective_all','FR_pref_all','Channel_depth_all','CellType_all','PSTH_all',...
                    'Hemisphere_all','Mice_all','t');
    Animals_list{I_animal}

    Psth_nostim_all={};
    Psth_stim_left_all={};
    Psth_stim_right_all={};
    Psth_nostim_mean_all={};
    Psth_stim_left_mean_all={};
    Psth_stim_right_mean_all={};
    Psth_useful_unit_hemisphere_all={};
    
    Psth_stimL_yes_all={};
    Psth_stimL_yes_mean_all={};
    Psth_stimL_no_all={};
    Psth_stimL_no_mean_all={};
    Psth_stimL_correct_all={};
    Psth_stimL_correct_mean_all={};
    Psth_stimL_error_all={};
    Psth_stimL_error_mean_all={};
    Psth_stimR_yes_all={};
    Psth_stimR_yes_mean_all={};
    Psth_stimR_no_all={};
    Psth_stimR_no_mean_all={};
    Psth_stimR_correct_all={};
    Psth_stimR_correct_mean_all={};
    Psth_stimR_error_all={};
    Psth_stimR_error_mean_all={};

    sel_nonstm_all={};
    sel_stim_leftALM_all={};
    sel_stim_rightALM_all={};
    sel_stim_Bilat_all={};
    sig_selective_Selected_all={};
    sel_unitlist_all={};
    sel_unit_hemisphere_all={};
   sig_selective_Selected_all2={};
    sel_unitlist_all2={};
    sel_unit_hemisphere_all2={};

    % 变量名含义为刺激侧
    sel_stimulus_left_all={};
    sel_stimulus_right_all={};
    sel_choice_left_all={}; %刺激左侧，记录右侧
    sel_choice_right_all={}; %刺激右侧，记录左侧
    sel_stimulus_nostim_all={};
    sel_choice_nostim_all={};

    % Nd变量名含义均是记录侧
    NdselLALM_all={};
    NdselRALM_all={};
    NdselLALM_stimulus_all={}; 
    NdselRALM_stimulus_all={}; 
    NdselLALM_choice_all={};
    NdselRALM_choice_all={};
    NsdrselLALM_all={};
    NsdrselRALM_all={};
    Mice_session_all=[];

    n_session=0;

    for i_session=1:length(obj.PSTH_all)
    
        i_session
        n_session=n_session+1;
    
        N_trials=obj.N_trials_all{i_session};
        sig_selective=obj.Sig_selective_all{i_session};         %[sample delay resposne sdr]
        FR_pref=obj.FR_pref_all{i_session};
        Channel_depth=obj.Channel_depth_all{i_session};
        Celltype=obj.CellType_all{i_session};
        PSTH=obj.PSTH_all{i_session};
        Hemisphere=obj.Hemisphere_all{i_session};  %1是左半球，2是右半球，确认于PSTH code
        Mouse_session=obj.Mice_all(i_session,:);
        t=obj.t;
    
        for j=1:size(PSTH,1)
            for k=1:size(PSTH,2)
                if isnan(mean(mean(PSTH{j,k})))
                PSTH{j,k}=[];
                end
            end
        end

        %opto tagging, identify putative inhibitory neurons strongly excited by light
        tsti=find(t>-1.7&t<-0.9);%early delay
        optotag_cell=zeros(size(PSTH,1),1);
        for j=1:size(PSTH,1)

            psthnostim=[PSTH{j,1};PSTH{j,2};PSTH{j,9};PSTH{j,10}];
            psthleft=[PSTH{j,3};PSTH{j,4};PSTH{j,11};PSTH{j,12}]; %刺激左
            psthright=[PSTH{j,5};PSTH{j,6};PSTH{j,13};PSTH{j,14}]; %刺激右
            psthbilat=[PSTH{j,7};PSTH{j,8};PSTH{j,15};PSTH{j,16}]; %刺激双侧

            if Hemisphere(j)==1&prod(size(psthleft))>0&prod(size(psthnostim))>0
                if mean(mean(psthleft(:,tsti)))>2*mean(mean(psthnostim(:,tsti)))|mean(mean(psthbilat(:,tsti)))>2*mean(mean(psthnostim(:,tsti)))
                    optotag_cell(j)=1;%strongly excited by light, putative inhibitory neuron
                end
            end
            if Hemisphere(j)==2&prod(size(psthright))>0&prod(size(psthnostim))>0
                if mean(mean(psthright(:,tsti)))>2*mean(mean(psthnostim(:,tsti)))|mean(mean(psthbilat(:,tsti)))>2*mean(mean(psthnostim(:,tsti)))
                    optotag_cell(j)=1;
                end
            end
        end

        sig_selective(isnan(sig_selective))=0;
        i_selective = (sig_selective(:,1)|sig_selective(:,2)|sig_selective(:,3)) & Celltype'==1 & (abs(FR_pref(:,1))>0.5|abs(FR_pref(:,2))>0.5|abs(FR_pref(:,3))>0.5)&optotag_cell==0;
        i_n_trial = sum(N_trials(1:size(PSTH,1),:),2)>100 & N_trials(1:size(PSTH,1),1)>5 & N_trials(1:size(PSTH,1),2)>5 & sum(N_trials(1:size(PSTH,1),[3 11]),2)>=2 & sum(N_trials(1:size(PSTH,1),[4 12]),2)>=2 & sum(N_trials(1:size(PSTH,1),[5 13]),2)>=2 & sum(N_trials(1:size(PSTH,1),[6 14]),2)>=2; 
        n_useful=find(Celltype'==1&i_n_trial&optotag_cell==0)';%only use excitatory neurons

        % ===== 计算  trialtype  selectivity =====
        sel_nonstm = [];
        sel_stim_leftALM = [];
        sel_stim_rightALM = [];
        sel_stim_Bilat = [];
        sig_selective_Selected = [];
        sel_unitlist=[];
        sel_unit_hemisphere=[];
        
        for i_cell =n_useful %use loose selective excitatory neurons
            
            if ~isempty([PSTH{i_cell,1};PSTH{i_cell,9}])&~isempty([PSTH{i_cell,2};PSTH{i_cell,10}])&...
                    ~isempty([PSTH{i_cell,3};PSTH{i_cell,11}])&~isempty([PSTH{i_cell,4};PSTH{i_cell,12}])&...
                    ~isempty([PSTH{i_cell,5};PSTH{i_cell,13}])&~isempty([PSTH{i_cell,6};PSTH{i_cell,14}])&...
                    size(PSTH{i_cell,9},1)>=2&size(PSTH{i_cell,10},1)>=2&...
                    size(PSTH{i_cell,11},1)>=2&size(PSTH{i_cell,12},1)>=2&...
                    size(PSTH{i_cell,13},1)>=2&size(PSTH{i_cell,14},1)>=2
                    
                psth_temp1=[PSTH{i_cell,1};PSTH{i_cell,9}];%lick right correct and error trials
                psth_temp2=[PSTH{i_cell,2};PSTH{i_cell,10}];%lick left correct and error trials
                correct_trial_n1=size(PSTH{i_cell,1},1);
                correct_trial_n2=size(PSTH{i_cell,2},1);
                
                FR_pref_screen=[];
                Half_sel_nonstm=[];
                
                for rep=1:30
                trialn1=randsample(size(psth_temp1,1),size(psth_temp1,1),0);
                trialn2=randsample(size(psth_temp2,1),size(psth_temp2,1),0);
                       
                trial_half11=trialn1(1:floor(length(trialn1)/2));
                trial_half12=trialn1(floor(length(trialn1)/2)+1:end);
                trial_half21=trialn2(1:floor(length(trialn2)/2));
                trial_half22=trialn2(floor(length(trialn2)/2)+1:end);      
                
                trial_label11=trial_half11(find(trial_half11<=correct_trial_n1));%screenCtrhalf only use correct trial as label
                trial_label12=trial_half12(find(trial_half12<=correct_trial_n1));
                trial_label21=trial_half21(find(trial_half21<=correct_trial_n2));
                trial_label22=trial_half22(find(trial_half22<=correct_trial_n2));
                
                FR_pref_screen1=mean(psth_temp1(trial_label11,:))-mean(psth_temp2(trial_label21,:));
                FR_pref_screen(end+1,:)=FR_pref_screen1;
                FR_pref_screen1=[mean(FR_pref_screen1(find(t<-1.7&t>-3))) mean(FR_pref_screen1(find(t<0&t>-1.7))) mean(FR_pref_screen1(find(t<1.3&t>0)))];
                W1=sign(mean(FR_pref_screen1));
                FR_pref_screen2=mean(psth_temp1(trial_label12,:))-mean(psth_temp2(trial_label22,:));
                FR_pref_screen(end+1,:)=FR_pref_screen2;
                FR_pref_screen2=[mean(FR_pref_screen2(find(t<-1.7&t>-3))) mean(FR_pref_screen2(find(t<0&t>-1.7))) mean(FR_pref_screen2(find(t<1.3&t>0)))];
                W2=sign(mean(FR_pref_screen2));
                
                Half_sel_nonstm(end+1,:)     = mean([W2.*(mean(psth_temp1(trial_half11,:))-mean(psth_temp2(trial_half21,:)));...
                                            W1.*(mean(psth_temp1(trial_half12,:))-mean(psth_temp2(trial_half22,:)))]);
                end
               
                FR_pref_screen=mean(FR_pref_screen);
                FR_pref_screen=[mean(FR_pref_screen(find(t<-1.7&t>-3))) mean(FR_pref_screen(find(t<0&t>-1.7))) mean(FR_pref_screen(find(t<1.3&t>0)))];
                W3=sign(mean(FR_pref_screen));
                       
                sel_stim_leftALM(end+1,:)   = W3*(mean([PSTH{i_cell,3};PSTH{i_cell,11}])-mean([PSTH{i_cell,4};PSTH{i_cell,12}]));
                sel_stim_rightALM(end+1,:)   = W3*(mean([PSTH{i_cell,5};PSTH{i_cell,13}])-mean([PSTH{i_cell,6};PSTH{i_cell,14}]));
                sel_stim_Bilat(end+1,:)      = W3*(mean([PSTH{i_cell,7};PSTH{i_cell,15}])-mean([PSTH{i_cell,8};PSTH{i_cell,16}]));
                sel_nonstm(end+1,:)     = mean(Half_sel_nonstm);
                
                sig_selective_Selected(end+1,:) = sig_selective(i_cell,:);
                sel_unitlist(end+1,:)=i_cell;
                sel_unit_hemisphere(end+1,:)=Hemisphere(i_cell);
        end
        end

        sel_nonstm_all{n_session}=sel_nonstm;
        sel_stim_leftALM_all{n_session}=sel_stim_leftALM;
        sel_stim_rightALM_all{n_session}=sel_stim_rightALM;
        sel_stim_Bilat_all{n_session}=sel_stim_Bilat;
        sig_selective_Selected_all{n_session}=sig_selective_Selected;
        sel_unitlist_all{n_session}=sel_unitlist;
        sel_unit_hemisphere_all{n_session}=sel_unit_hemisphere;
        
        %find delay selective neurons
        NdselLALM=find(sig_selective(sel_unitlist,2)==1&abs(FR_pref(sel_unitlist,2))>0.5&...
            sum(N_trials(sel_unitlist,[5 13]),2)>=5&sum(N_trials(sel_unitlist,[6 14]),2)>=5&sel_unit_hemisphere==1);%delay selective
        NdselRALM=find(sig_selective(sel_unitlist,2)==1&abs(FR_pref(sel_unitlist,2))>0.5&...
            sum(N_trials(sel_unitlist,[3 11]),2)>=5&sum(N_trials(sel_unitlist,[4 12]),2)>=5&sel_unit_hemisphere==2);%delay selective
        NdselLALM_all{n_session}=NdselLALM;
        NdselRALM_all{n_session}=NdselRALM;
        
        %find sample/delay/response epoch selective neurons
        tmp=find(abs(FR_pref(sel_unitlist,1))>0.5|abs(FR_pref(sel_unitlist,2))>0.5|abs(FR_pref(sel_unitlist,3))>0.5);%sample/delay/response selective
        NsdrselLALM=intersect(tmp,find(sel_unit_hemisphere==1));
        NsdrselRALM=intersect(tmp,find(sel_unit_hemisphere==2));
        NsdrselLALM_all{n_session}=NsdrselLALM;
        NsdrselRALM_all{n_session}=NsdrselRALM;
        Mice_session_all(n_session,:)=Mouse_session;


sel_stimulus_nostim = [];
sel_stimulus_left = [];
sel_stimulus_right = [];
sel_choice_nostim = [];
sel_choice_left = [];
sel_choice_right = [];
sig_selective_Selected2 = [];
sel_unitlist2=[];
sel_unit_hemisphere2=[];
for i_cell =intersect(n_useful,find(i_selective & i_n_trial)')   
    if ~isempty([PSTH{i_cell,1};PSTH{i_cell,9}])&~isempty([PSTH{i_cell,2};PSTH{i_cell,10}])%&...
        if size(PSTH{i_cell,9},1)>=2&size(PSTH{i_cell,10},1)>=2
            
            %for stimulus selectivity
            stimulus_w_time_tmp=[];%prefered direction, w>0 prefer right, w<0 prefer left
            stimulus_sel_nonstm_tmp = [];%selectivity
            for i_time=1:51%bin step 100ms, 51 bins
                tbin=100;
                psth1=[];
                psth2=[];
                psth3=[];
                psth4=[];
                
                if ~isempty(PSTH{i_cell,1})
                    psth1=mean(PSTH{i_cell,1}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,10})
                    psth2=mean(PSTH{i_cell,10}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,9})
                    psth3=mean(PSTH{i_cell,9}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,2})
                    psth4=mean(PSTH{i_cell,2}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                w=sign(mean(psth1)-mean(psth2)+mean(psth3)-mean(psth4));
                stimulus_w_time_tmp(i_time)=w;
            end

             FR_pref_screen = [];
             Half_sel_nonstm = [];
            for ktmp=1:10
                krand=randsample(size(PSTH{i_cell,1},1),size(PSTH{i_cell,1},1),0);
                krand11=krand(1:floor(length(krand)/2));
                krand12=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,10},1),size(PSTH{i_cell,10},1),0);
                krand21=krand(1:floor(length(krand)/2));
                krand22=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,9},1),size(PSTH{i_cell,9},1),0);
                krand31=krand(1:floor(length(krand)/2));
                krand32=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,2},1),size(PSTH{i_cell,2},1),0);
                krand41=krand(1:floor(length(krand)/2));
                krand42=krand(floor(length(krand)/2)+1:end);

                w1tmp=mean(PSTH{i_cell,1}(krand11,:),1)-mean(PSTH{i_cell,10}(krand21,:),1)+mean(PSTH{i_cell,9}(krand31,:),1)-mean(PSTH{i_cell,2}(krand41,:),1);
                FR_pref_screen(end+1,:)=w1tmp;
                w1tmp=[mean(w1tmp(find(t>-3&t<-1.7))) mean(w1tmp(find(t>-1.7&t<0))) mean(w1tmp(find(t>0&t<1.5)))];
                w1=sign(mean(w1tmp));

                w2tmp=mean(PSTH{i_cell,1}(krand12,:),1)-mean(PSTH{i_cell,10}(krand22,:),1)+mean(PSTH{i_cell,9}(krand32,:),1)-mean(PSTH{i_cell,2}(krand42,:),1);
                FR_pref_screen(end+1,:)=w2tmp;
                w2tmp=[mean(w2tmp(find(t>-3&t<-1.7))) mean(w2tmp(find(t>-1.7&t<0))) mean(w2tmp(find(t>0&t<1.5)))];
                w2=sign(mean(w2tmp));

                sel_half = [w2*(mean(PSTH{i_cell,1}(krand11,:),1)-mean(PSTH{i_cell,10}(krand21,:),1)+...
                                mean(PSTH{i_cell,9}(krand31,:),1)-mean(PSTH{i_cell,2}(krand41,:),1))/2;
                                w1*(mean(PSTH{i_cell,1}(krand12,:),1)-mean(PSTH{i_cell,10}(krand22,:),1)+...
                                mean(PSTH{i_cell,9}(krand32,:),1)-mean(PSTH{i_cell,2}(krand42,:),1))/2];
                Half_sel_nonstm(end+1,:) = mean(sel_half,1); 
                stimulus_sel_nonstm_tmp = [stimulus_sel_nonstm_tmp; sel_half];
            end
            FR_pref_screen=mean(FR_pref_screen);
            FR_pref_screen=[nanmean(FR_pref_screen(find(t<-1.7&t>-3))) nanmean(FR_pref_screen(find(t<0&t>-1.7))) nanmean(FR_pref_screen(find(t<1.3&t>0)))];
            W_stimulus = sign(nanmean(FR_pref_screen));     
            sel_stimulus_nostim(end+1,:) = mean(Half_sel_nonstm);    %stimulus selectivity
            sel_stimulus_left(end+1,:)   = W_stimulus* ((mean(PSTH{i_cell,3}) + mean(PSTH{i_cell,11}) - mean(PSTH{i_cell,4}) - mean(PSTH{i_cell,12}) )/2);
            sel_stimulus_right(end+1,:)   = W_stimulus* ((mean(PSTH{i_cell,5}) + mean(PSTH{i_cell,13}) - mean(PSTH{i_cell,6}) - mean(PSTH{i_cell,14}) )/2);

            %for choice selectivity
            choice_w_time_tmp=[];
            choice_sel_nonstm_tmp = [];
            for i_time=1:51%timewindow 100ms 51 bins 
                tbin=100;
                psth1=[];
                psth2=[];
                psth3=[];
                psth4=[];
                if ~isempty(PSTH{i_cell,1})
                    psth1=mean(PSTH{i_cell,1}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,9})
                    psth2=mean(PSTH{i_cell,9}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,10})
                    psth3=mean(PSTH{i_cell,10}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,2})
                    psth4=mean(PSTH{i_cell,2}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                w=sign(mean(psth1)-mean(psth2)+mean(psth3)-mean(psth4));
                choice_w_time_tmp(i_time)=w;
            end

            FR_pref_screen = [];
            Half_sel_nonstm = [];
            for ktmp=1:10
                krand=randsample(size(PSTH{i_cell,1},1),size(PSTH{i_cell,1},1),0);
                krand11=krand(1:floor(length(krand)/2));
                krand12=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,9},1),size(PSTH{i_cell,9},1),0);
                krand21=krand(1:floor(length(krand)/2));
                krand22=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,10},1),size(PSTH{i_cell,10},1),0);
                krand31=krand(1:floor(length(krand)/2));
                krand32=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,2},1),size(PSTH{i_cell,2},1),0);
                krand41=krand(1:floor(length(krand)/2));
                krand42=krand(floor(length(krand)/2)+1:end);
                w1tmp=mean(PSTH{i_cell,1}(krand11,:),1)-mean(PSTH{i_cell,9}(krand21,:),1)+mean(PSTH{i_cell,10}(krand31,:),1)-mean(PSTH{i_cell,2}(krand41,:),1);
                FR_pref_screen(end+1,:)=w1tmp;
                w1tmp=[mean(w1tmp(find(t>-3&t<-1.7))) mean(w1tmp(find(t>-1.7&t<0))) mean(w1tmp(find(t>0&t<1.5)))];
                w1=sign(mean(w1tmp));
                
                w2tmp=mean(PSTH{i_cell,1}(krand12,:),1)-mean(PSTH{i_cell,9}(krand22,:),1)+mean(PSTH{i_cell,10}(krand32,:),1)-mean(PSTH{i_cell,2}(krand42,:),1);
                FR_pref_screen(end+1,:)=w2tmp;
                w2tmp=[mean(w2tmp(find(t>-3&t<-1.7))) mean(w2tmp(find(t>-1.7&t<0))) mean(w2tmp(find(t>0&t<1.5)))];
                w2=sign(mean(w2tmp));

                sel_half = [w2*(mean(PSTH{i_cell,1}(krand11,:),1)-mean(PSTH{i_cell,9}(krand21,:),1)+...
                                mean(PSTH{i_cell,10}(krand31,:),1)-mean(PSTH{i_cell,2}(krand41,:),1))/2;
                                w1*(mean(PSTH{i_cell,1}(krand12,:),1)-mean(PSTH{i_cell,9}(krand22,:),1)+...
                                mean(PSTH{i_cell,10}(krand32,:),1)-mean(PSTH{i_cell,2}(krand42,:),1))/2];
                Half_sel_nonstm(end+1,:) = mean(sel_half,1); 
                choice_sel_nonstm_tmp = [choice_sel_nonstm_tmp; sel_half];

            end
            FR_pref_screen=mean(FR_pref_screen);   
            FR_pref_screen = [ mean(FR_pref_screen(t>-3   & t<-1.7)), mean(FR_pref_screen(t>-1.7 & t<0)), mean(FR_pref_screen(t>0& t<1.5)) ];
            W_choice = sign(nanmean(FR_pref_screen));       
            sel_choice_nostim(end+1,:) = mean(Half_sel_nonstm);  %choice selectivity
            sel_choice_left(end+1,:)   = W_choice* ((mean(PSTH{i_cell,3}) + mean(PSTH{i_cell,12}) - mean(PSTH{i_cell,4}) - mean(PSTH{i_cell,11}) )/2);
            sel_choice_right(end+1,:)   = W_choice* ((mean(PSTH{i_cell,5}) + mean(PSTH{i_cell,14}) - mean(PSTH{i_cell,6}) - mean(PSTH{i_cell,13}) )/2);

        sig_selective_Selected2(end+1,:) = sig_selective(i_cell,:);
        sel_unitlist2(end+1,:)           = i_cell;
        sel_unit_hemisphere2(end+1,:)    = Hemisphere(i_cell);
        end
    end
end

        sig_selective_Selected_all2{n_session}=sig_selective_Selected2;
        sel_unitlist_all2{n_session}=sel_unitlist2;
        sel_unit_hemisphere_all2{n_session}=sel_unit_hemisphere2;

        sel_stimulus_left_all{n_session}=sel_stimulus_left;
        sel_stimulus_right_all{n_session}=sel_stimulus_right;
        sel_stimulus_nostim_all{n_session}=sel_stimulus_nostim;
        sel_choice_left_all{n_session}=sel_choice_left;
        sel_choice_right_all{n_session}=sel_choice_right;
        sel_choice_nostim_all{n_session}=sel_choice_nostim;

        %find delay selective neurons 
        NdselLALM_stimulus=find(sig_selective(sel_unitlist2,2)==1&abs(FR_pref(sel_unitlist2,2))>0.5&...
                N_trials(sel_unitlist2,5) >= 2 & N_trials(sel_unitlist2,13) >= 2 & ...
                N_trials(sel_unitlist2,6) >= 2 & N_trials(sel_unitlist2,14) >= 2 & ...
                sel_unit_hemisphere2==1);%左半球神经元
        NdselRALM_stimulus=find(sig_selective(sel_unitlist2,2)==1&abs(FR_pref(sel_unitlist2,2))>0.5&...
                N_trials(sel_unitlist2,3) >= 2 & N_trials(sel_unitlist2,11) >= 2 & ...
                N_trials(sel_unitlist2,4) >= 2 & N_trials(sel_unitlist2,12) >= 2 & ...
                sel_unit_hemisphere2==2);

        NdselLALM_choice=find(sig_selective(sel_unitlist2,2)==1&abs(FR_pref(sel_unitlist2,2))>0.5&...
                N_trials(sel_unitlist2,5) >= 2 & N_trials(sel_unitlist2,13) >= 2 & ...
                N_trials(sel_unitlist2,6) >= 2 & N_trials(sel_unitlist2,14) >= 2 & ...
                sel_unit_hemisphere2==1);
        NdselRALM_choice=find(sig_selective(sel_unitlist2,2)==1&abs(FR_pref(sel_unitlist2,2))>0.5&...
                N_trials(sel_unitlist2,3) >= 2 & N_trials(sel_unitlist2,11) >= 2 & ...
                N_trials(sel_unitlist2,4) >= 2 & N_trials(sel_unitlist2,12) >= 2 & ...
                sel_unit_hemisphere2==2);

        NdselLALM_stimulus_all{n_session}=NdselLALM_stimulus; %左侧ALM的记录
        NdselRALM_stimulus_all{n_session}=NdselRALM_stimulus;
        NdselLALM_choice_all{n_session}=NdselLALM_choice;
        NdselRALM_choice_all{n_session}=NdselRALM_choice;

        %find sample/delay/response epoch selective neurons
        tmp=find(abs(FR_pref(sel_unitlist,1))>0.5|abs(FR_pref(sel_unitlist,2))>0.5|abs(FR_pref(sel_unitlist,3))>0.5);%sample/delay/response selective
        NsdrselLALM=intersect(tmp,find(sel_unit_hemisphere==1));
        NsdrselRALM=intersect(tmp,find(sel_unit_hemisphere==2));
        NsdrselLALM_all{n_session}=NsdrselLALM;
        NsdrselRALM_all{n_session}=NsdrselRALM;
        Mice_session_all(n_session,:)=Mouse_session;
    end

    clear obj
    save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_',Animals_list{I_animal},'.mat'],...
            'sel_nonstm_all','sel_stim_leftALM_all','sel_stim_rightALM_all','sel_stim_Bilat_all',...
            'sel_stimulus_left_all','sel_stimulus_right_all','sel_choice_right_all', 'sel_choice_left_all','sel_stimulus_nostim_all','sel_choice_nostim_all',...
            'sig_selective_Selected_all','sel_unitlist_all','sel_unit_hemisphere_all',...
            'NdselLALM_stimulus_all','NdselRALM_stimulus_all', 'NdselLALM_choice_all','NdselRALM_choice_all', 'NdselLALM_all','NdselRALM_all', ...
            'NsdrselLALM_all','NsdrselRALM_all', 'Mice_session_all','t','tsti');
    clearvars -except Animals_list I_animal N_animals filepath
end

%% combine all sessions selectivity related variables and save them to one file
clear all
close all
addpath('.\helper_functions\')
filepath='C:\Users\XinHao\Desktop';

% list of animals to analyze
Animals_list = {
    % standard task
    'BAYLORGC4';...  % VGAT
    'BAYLORGC12';...  % PVReachR
    'BAYLORGC13';...  % PVReachR
    'BAYLORGC15';...  % PVReachR
    'BAYLORGC17';...  % PVReachR    
    'BAYLORGC18';...  % PVReachR
    'BAYLORGC19';...  % PVReachR
    'BAYLORGC21';...  % PVReachR
    'BAYLORGC22';...  % PVReachR
    'BAYLORGC25';...  % PVReachR
    'BAYLORGC26';...  % PVReachR
    'BAYLORGC50';...  % VGAT
    'BAYLORGC72';...  % VGAT
    'BAYLORGC75';...  % VGAT
    'BAYLORGC86';...  % PVReachR
    'BAYLORGC87';...  % PVReachR
    'BAYLORGC93';...  % PVReachR
    'BAYLORGC95';...  % PVReachR
    
    % reversed contingency task
    'BAYLORGC61';...  % VGAT
    'BAYLORGC69';...  % VGAT
    'BAYLORGC71';...  % VGAT
    'BAYLORGC74';...  % VGAT
    'BAYLORGC78';...  % VGAT
    
    % fully reversed task
    'BAYLORGC79';...  % PVReachR
    'BAYLORGC80';...  % PVReachR
    'BAYLORGC84';...  % PVReachR
    'BAYLORGC85';...  % PVReachR
    'BAYLORGC83';...  % PVReachR
    'BAYLORGC88';...  % PVReachR
    'BAYLORGC81';...  % PVReachR
    'BAYLORGC89';...  % PVReachR
    'BAYLORGC90';...  % PVReachR
    'BAYLORGC94';...  % PVReachR
    
    % reversed tactile stimulus task
    'BAYLORGC91';...  % PVReachR
    'BAYLORGC92';...  % PVReachR
    'BAYLORGC96';...  % PVReachR
    'BAYLORGC98';...  % PVReachR
    'BAYLORGC99';...  % PVReachR
    'BAYLORGC100';...  % PVReachR  
    };

Psth_nostim_all={};
Psth_stim_left_all={};
Psth_stim_right_all={};
Psth_nostim_mean_all={};
Psth_stim_left_mean_all={};
Psth_stim_right_mean_all={};
Psth_useful_unit_hemisphere_all={};

sel_nonstm_all={};
sel_stim_leftALM_all={};
sel_stim_rightALM_all={};
sel_stim_Bilat_all={};
sig_selective_Selected_all={};
sel_unitlist_all={};
sel_unit_hemisphere_all={};
sel_stimulus_left_all={};
sel_stimulus_right_all={};
sel_choice_left_all={};
sel_choice_right_all={};
sel_stimulus_nostim_all={};
sel_choice_nostim_all={};

NsdrselLALM_all={};
NsdrselRALM_all={};
Mice_session_all=[];
NdselLALM_stimulus_all={};
NdselRALM_stimulus_all={};
NdselLALM_choice_all={};
NdselRALM_choice_all={};

n_session=0;
N_animals = length(Animals_list); 

for I_animal = 1:N_animals
    Animals_list{I_animal};
    obj=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_',Animals_list{I_animal},'.mat']);

    try
        for i_session=1:length(obj.sel_nonstm_all)
            i_session;
            n_session=n_session+1;
    
            sel_stimulus_left_all{n_session} = obj.sel_stimulus_left_all{i_session};
            sel_stimulus_right_all{n_session}= obj.sel_stimulus_right_all{i_session};
            sel_stimulus_nostim_all{n_session}= obj.sel_stimulus_nostim_all{i_session};
            sel_choice_left_all{n_session}= obj.sel_choice_left_all{i_session};
            sel_choice_right_all{n_session}= obj.sel_choice_right_all{i_session};
            sel_choice_nostim_all{n_session}= obj.sel_choice_nostim_all{i_session};

            sel_nonstm_all{n_session}=obj.sel_nonstm_all{i_session};
            sel_stim_leftALM_all{n_session}=obj.sel_stim_leftALM_all{i_session};
            sel_stim_rightALM_all{n_session}=obj.sel_stim_rightALM_all{i_session};
            sel_stim_Bilat_all{n_session}=obj.sel_stim_Bilat_all{i_session};
            sig_selective_Selected_all{n_session}=obj.sig_selective_Selected_all{i_session};
            sel_unitlist_all{n_session}=obj.sel_unitlist_all{i_session};
            sel_unit_hemisphere_all{n_session}=obj.sel_unit_hemisphere_all{i_session};
    
            NdselLALM_all{n_session}=obj.NdselLALM_all{i_session};
            NdselRALM_all{n_session}=obj.NdselRALM_all{i_session};
            NdselLALM_stimulus_all{n_session}=obj.NdselLALM_stimulus_all{i_session};
            NdselRALM_stimulus_all{n_session}=obj.NdselRALM_stimulus_all{i_session};
            NdselLALM_choice_all{n_session}=obj.NdselLALM_choice_all{i_session};
            NdselRALM_choice_all{n_session}=obj.NdselRALM_choice_all{i_session};

            NsdrselLALM_all{n_session}=obj.NsdrselLALM_all{i_session};
            NsdrselRALM_all{n_session}=obj.NsdrselRALM_all{i_session};
            Mice_session_all(n_session,:)=obj.Mice_session_all(i_session,:);
        end

    catch ME
        warning('Error in session %d of %s: %s', i_session, Animals_list{I_animal}, ME.message);
    end
end

t=obj.t;
tsti=obj.tsti;
clear obj
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat'],...
        'sel_nonstm_all','sel_stim_leftALM_all','sel_stim_rightALM_all','sel_stim_Bilat_all',...
        'sel_stimulus_left_all', 'sel_stimulus_right_all', 'sel_choice_left_all', 'sel_choice_right_all', 'sel_stimulus_nostim_all', 'sel_choice_nostim_all',...
        'sig_selective_Selected_all','sel_unitlist_all','sel_unit_hemisphere_all',...
        'NdselLALM_stimulus_all','NdselRALM_stimulus_all','NdselLALM_choice_all','NdselRALM_choice_all', 'NdselLALM_all','NdselRALM_all',...
        'NsdrselLALM_all','NsdrselRALM_all', 'Mice_session_all','t','tsti');

%%  原文章的图，不区分stimulus and choice
clear;
filepath='C:\Users\XinHao\Desktop';
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
figure;
for m=1:2
    if m==1
       N=47; %modular session 47 BAYLORGC86_2019_07_31; asym session 35 BAYLORGC50_2018_12_07
    else
        N=35;
    end

    for j=1:2
        sel_nonstm=[];
        sel_sti=[];
        switch j
            case 1
                sel_nonstm=[sel_nonstm;sel_nonstm_all{N}(NdselLALM_all{N},:)];
                sel_sti=[sel_sti;sel_stim_rightALM_all{N}(NdselLALM_all{N},:)];
            case 2
                sel_nonstm=[sel_nonstm;sel_nonstm_all{N}(NdselRALM_all{N},:)];
                sel_sti=[sel_sti;sel_stim_leftALM_all{N}(NdselRALM_all{N},:)];
        end
        subplot(3,3,(m-1)*3+j);
        [hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
        outlinebounds(hl,hp);
        hold on;
        [hl,hp]=boundedline(t,mean(sel_sti),std(sel_sti)/sqrt(size(sel_sti,1)),'-b','transparency',0.1);
        outlinebounds(hl,hp);
        hold on;
        Ylim=round(1.3*max([max(mean(sel_sti)) max(mean(sel_nonstm))]));
        ylim([-1, 4]);
        line([0 0],[-1 Ylim],'color','k')
        line([-1.7 -1.7],[-1 Ylim],'color','k')
        line([-3 -3],[-1 Ylim],'color','k')
        line([min(t) max(t)],[0 0],'color','k')
        xlabel('time (s)')
        if j==1
            ylabel('LALM selectivity(spikes/s)');
        else
            ylabel('RALM selectivity(spikes/s)');
        end   

        Nsel=size(sel_nonstm,1);
        if j==1
            title(['StiRight',' Nsel=',num2str(Nsel)]);
        else
            title(['StiLeft',' Nsel=',num2str(Nsel)]);
        end
        Si=[];
        for k=1:1000
            krand=randsample(size(sel_sti,1),size(sel_sti,1),1);
            Si(k)=mean(mean(sel_sti(krand,find(t>-1.7&t<-0.9))))/mean(mean(sel_nonstm(krand,find(t>-1.7&t<-0.9))));
        end
        subplot(3,3,(m-1)*3+3);hold on;
        bar(j,mean(Si));hold on;
        errorbar(j,mean(Si),std(Si));hold on;
        ylabel('modularity');
        xlabel('Left Right (ALM)');
    end
end

% 图 7: LALM  平均
subplot(3,3,7);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_nonstm_all)
    if ~isempty(NdselLALM_stimulus_all{N})
        sel_nonstm = [sel_nonstm; sel_nonstm_all{N}(NdselLALM_all{N}, :)];
        sel_sti    = [sel_sti; sel_stim_rightALM_all{N}(NdselLALM_all{N}, :)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), ...
    std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), ...
    std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);
line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('LALM all neurons sel (spikes/s)');
ylim([-1 4]);

% 图 8: RALM neuron-level 平均
subplot(3,3,8);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_nonstm_all)
    if ~isempty(NdselRALM_stimulus_all{N})
        sel_nonstm = [sel_nonstm; sel_nonstm_all{N}(NdselRALM_all{N}, :)];
        sel_sti    = [sel_sti; sel_stim_leftALM_all{N}(NdselRALM_all{N}, :)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), ...
    std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), ...
    std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);
line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('RALM all neurons sel (spikes/s)');
ylim([-1 4]);

saveas(gcf,'FigureTrialtype.fig','fig');

%% stimulus图
clear;
filepath='C:\Users\XinHao\Desktop';
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
figure;
 for m=1:2
     if m==1
        N=47; %modular session 47 BAYLORGC86_2019_07_31; asym session 35 BAYLORGC50_2018_12_07
     else
         N=35;
     end
     for j=1:2
         sel_nonstm=[];
         sel_sti=[];
         switch j
             case 1
                 sel_nonstm=[sel_nonstm;sel_stimulus_nostim_all{N}(NdselLALM_stimulus_all{N},:)]; %无刺激条件下的左侧神经元
                 sel_sti=[sel_sti;sel_stimulus_right_all{N}(NdselLALM_stimulus_all{N},:)]; %右侧刺激条件下的左侧神经元
             case 2
                 sel_nonstm=[sel_nonstm;sel_stimulus_nostim_all{N}(NdselRALM_stimulus_all{N},:)];
                 sel_sti=[sel_sti;sel_stimulus_left_all{N}(NdselRALM_stimulus_all{N},:)];
         end
         subplot(3,3,(m-1)*3+j);
         [hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
         outlinebounds(hl,hp);
         hold on;
         [hl,hp]=boundedline(t,mean(sel_sti),std(sel_sti)/sqrt(size(sel_sti,1)),'-b','transparency',0.1);
         outlinebounds(hl,hp);
         hold on;
         Ylim=round(1.3*max([max(mean(sel_sti)) max(mean(sel_nonstm))]));
         ylim([-1, 4]);
         line([0 0],[-1 Ylim],'color','k')
         line([-1.7 -1.7],[-1 Ylim],'color','k')
         line([-3 -3],[-1 Ylim],'color','k')
         line([min(t) max(t)],[0 0],'color','k')
         xlabel('time (s)')
         if j==1
             ylabel('LALM selectivity(spikes/s)');
         else
             ylabel('RALM selectivity(spikes/s)');
         end   
 
         Nsel=size(sel_nonstm,1);
         if j==1
             title(['RALMsti',' Nsel=',num2str(Nsel)]);
         else
             title(['LALMsti',' Nsel=',num2str(Nsel)]);
         end
         Si=[];
         for k=1:1000
             krand=randsample(size(sel_sti,1),size(sel_sti,1),1);
             Si(k)=mean(mean(sel_sti(krand,find(t>-1.7&t<-0.9))))/mean(mean(sel_nonstm(krand,find(t>-1.7&t<-0.9))));
         end
         subplot(3,3,(m-1)*3+3);hold on;
         bar(j,mean(Si));hold on;
         errorbar(j,mean(Si),std(Si));hold on;
         ylabel('modularity');
         xlabel('Left Right (ALM)');
     end
 end

% 图 7: LALM 平均
subplot(3,3,7);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_stimulus_nostim_all)
    if ~isempty(NdselLALM_stimulus_all{N})
        sel_nonstm = [sel_nonstm; sel_stimulus_nostim_all{N}(NdselLALM_stimulus_all{N}, :)];
        sel_sti    = [sel_sti; sel_stimulus_right_all{N}(NdselLALM_stimulus_all{N}, :)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);

line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('LALM all neurons sel (spikes/s)');
ylim([-1 4]);


% 图 8: RALM neuron-level 平均（剔除含 NaN 的 neuron）
subplot(3,3,8);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_stimulus_nostim_all)
    if ~isempty(NdselRALM_stimulus_all{N})
        sel_nonstm = [sel_nonstm; sel_stimulus_nostim_all{N}(NdselRALM_stimulus_all{N}, :)];
        sel_sti    = [sel_sti; sel_stimulus_left_all{N}(NdselRALM_stimulus_all{N}, :)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);

line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('RALM all neurons sel (spikes/s)');
ylim([-1 4]);

saveas(gcf,'FigureStimulus.fig','fig');

%% choice figure
clear;
filepath='C:\Users\XinHao\Desktop';
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
figure;
set(gcf,'Color','w');
for m=1:2
    if m==1
       N=47; %47,35
    else
        N=35;
    end
for j=1:2
    sel_nonstm=[];
    sel_sti=[];
    switch j
        case 1
            sel_nonstm=[sel_nonstm;sel_choice_nostim_all{N}(NdselLALM_choice_all{N},:)];
            sel_sti=[sel_sti;sel_choice_right_all{N}(NdselLALM_choice_all{N},:)];
        case 2
            sel_nonstm=[sel_nonstm;sel_choice_nostim_all{N}(NdselRALM_choice_all{N},:)];
            sel_sti=[sel_sti;sel_choice_left_all{N}(NdselRALM_choice_all{N},:)];
    end
subplot(3,3,(m-1)*3+j);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_sti),std(sel_sti,[],1)/sqrt(size(sel_sti,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_sti)) max(mean(sel_nonstm))]));
ylim([-1, 4]);
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
if j==1
ylabel('LALM selectivity(spikes/s)');
else if j==2
ylabel('RALM selectivity(spikes/s)');
    end   
end
Nsel=size(sel_nonstm,1);
if j==1
    title(['RALMcho',' Nsel=',num2str(Nsel)],'Color','k');
end
if j==2
    title(['LALMcho',' Nsel=',num2str(Nsel)],'Color','k');
end
Si=[];
for k=1:1000
    krand=randsample(size(sel_sti,1),size(sel_sti,1),1);
    Si(k)=mean(mean(sel_sti(krand,find(t>-1.7&t<-0.9))))/mean(mean(sel_nonstm(krand,find(t>-1.7&t<-0.9))));
end
subplot(3,3,(m-1)*3+3);hold on;
bar(j,mean(Si));hold on;
errorbar(j,mean(Si),std(Si));hold on;
ylabel('modularity');
xlabel('Left Right (ALM)');
end
end

% 图 7: LALM neuron-level 平均（剔除含 NaN 的 neuron）
subplot(3,3,7);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_choice_nostim_all)
    if ~isempty(NdselLALM_choice_all{N})
        tmp_nonstm = sel_choice_nostim_all{N}(NdselLALM_choice_all{N}, :);
        tmp_sti    = sel_choice_right_all{N}(NdselLALM_choice_all{N}, :);
        valid_idx  = all(~isnan(tmp_nonstm),2) & all(~isnan(tmp_sti),2);
        sel_nonstm = [sel_nonstm; tmp_nonstm(valid_idx,:)];
        sel_sti    = [sel_sti; tmp_sti(valid_idx,:)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);

line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('LALM all neurons sel (spikes/s)');
ylim([-1 4]);


% 图 8: RALM neuron-level 平均（剔除含 NaN 的 neuron）
subplot(3,3,8);
sel_nonstm = [];
sel_sti = [];
for N = 1:length(sel_choice_nostim_all)
    if ~isempty(NdselRALM_choice_all{N})
        tmp_nonstm = sel_choice_nostim_all{N}(NdselRALM_choice_all{N}, :);
        tmp_sti    = sel_choice_left_all{N}(NdselRALM_choice_all{N}, :);
        valid_idx  = all(~isnan(tmp_nonstm),2) & all(~isnan(tmp_sti),2);
        sel_nonstm = [sel_nonstm; tmp_nonstm(valid_idx,:)];
        sel_sti    = [sel_sti; tmp_sti(valid_idx,:)];
    end
end
[hl,hp] = boundedline(t, mean(sel_nonstm,1), std(sel_nonstm,[],1)/sqrt(size(sel_nonstm,1)), '-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp] = boundedline(t, mean(sel_sti,1), std(sel_sti,[],1)/sqrt(size(sel_sti,1)), '-b','transparency',0.1);
outlinebounds(hl,hp);

line([0 0],[-1 4],'color','k');
line([-1.7 -1.7],[-1 4],'color','k');
line([-3 -3],[-1 4],'color','k');
line([min(t) max(t)], [0 0], 'color', 'k');
xlabel('time (s)');
ylabel('RALM all neurons sel (spikes/s)');
ylim([-1 4]);


saveas(gcf,'FigureChoice.fig','fig');

%% delay selectivity, modularity, robustness, and other meta information
clear;
filepath='C:\Users\XinHao\Desktop';

selectivity_meta_str={'mice','session','perf_ctr','perf_ctr_Rtrial','perf_ctr_Ltrial',...
    'perf_stiLeft','perf_stiRight','perf_stiBilat','num_dselLALM','num_dselRALM',...
    'LALM_earlyD_sel','LALM_lateD_sel','RALM_earlyD_sel','RALM_lateD_sel',...
    'LALM_modularity','RALM_modularity','LALMstiLSelrec','LALMstiRSelrec',...
    'RALMstiLSelrec','RALMstiRSelrec','LRALMstiLSelrec','LRALMstiRSelrec',...
    'craniotomy_quality','penetration_quality','task_type','strain','sex','age','trainingdays',...
    'LALM_rate_nonstm','LALM_rate_stim_left','LALM_rate_stim_right',...
    'RALM_rate_nonstm','RALM_rate_stim_left','RALM_rate_stim_right'};
twin1=[-1.7 -0.9];
twin2=[-0.8 0];

obj1 = load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_behavior.mat']);
obj2 = load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
t = obj2.t;

selectivity_meta = [];

for j = 1:size(obj2.Mice_session_all,1)
    j;
    selectivity_meta(j,1:2) = obj2.Mice_session_all(j,1:2); % mice session name

    % behavior performance
    if obj1.behavior_performance_allSession(j,1)==obj2.Mice_session_all(j,1) && ...
       obj1.behavior_performance_allSession(j,2)==obj2.Mice_session_all(j,2)
        selectivity_meta(j,3:8) = obj1.behavior_performance_allSession(j,3:end);
    end

    % 取各类变量
    NdselLALM = obj2.NdselLALM_all{j};
    NdselRALM = obj2.NdselRALM_all{j};
    sel_nonstm = obj2.sel_nonstm_all{j};
    sel_stim_leftALM = obj2.sel_stim_leftALM_all{j};
    sel_stim_rightALM = obj2.sel_stim_rightALM_all{j};
    sel_stimulus_left = obj2.sel_stimulus_left_all{j};
    sel_stimulus_right = obj2.sel_stimulus_right_all{j};
    sel_choice_left = obj2.sel_choice_left_all{j};
    sel_choice_right = obj2.sel_choice_right_all{j};
    sel_stimulus_nostim = obj2.sel_stimulus_nostim_all{j};
    sel_choice_nostim = obj2.sel_choice_nostim_all{j};
    NdselLALM_stimulus = obj2.NdselLALM_stimulus_all{j};
    NdselRALM_stimulus = obj2.NdselRALM_stimulus_all{j};
    NdselLALM_choice = obj2.NdselLALM_choice_all{j};
    NdselRALM_choice = obj2.NdselRALM_choice_all{j};

    % number of delay selective neurons
    if isempty(NdselLALM)
        selectivity_meta(j,9) = NaN;
    else
        selectivity_meta(j,9) = length(NdselLALM);
    end
    if isempty(NdselRALM)
        selectivity_meta(j,10) = NaN;
    else
        selectivity_meta(j,10) = length(NdselRALM);
    end

    % 定义时间区间索引
    cols1 = find(t>twin1(1) & t<twin1(2));
    cols2 = find(t>twin2(1) & t<twin2(2));

    % LALM early/late delay
    if ~isempty(NdselLALM) && ~isempty(sel_nonstm)
        selectivity_meta(j,11) = mean(sel_nonstm(NdselLALM,cols1),'all','omitnan');
        selectivity_meta(j,12) = mean(sel_nonstm(NdselLALM,cols2),'all','omitnan');
    else
        selectivity_meta(j,11:12) = NaN;
    end

    % RALM early/late delay
    if ~isempty(NdselRALM) && ~isempty(sel_nonstm)
        selectivity_meta(j,13) = mean(sel_nonstm(NdselRALM,cols1),'all','omitnan');
        selectivity_meta(j,14) = mean(sel_nonstm(NdselRALM,cols2),'all','omitnan');
    else
        selectivity_meta(j,13:14) = NaN;
    end

    % trialtype modularity
    if ~isempty(NdselLALM) && ~isempty(sel_stim_rightALM) && ~isempty(sel_nonstm)
        selectivity_meta(j,15) = mean(sel_stim_rightALM(NdselLALM,cols1),'all','omitnan') / ...
                                 mean(sel_nonstm(NdselLALM,cols1),'all','omitnan');
    else
        selectivity_meta(j,15) = NaN;
    end
    if ~isempty(NdselRALM) && ~isempty(sel_stim_leftALM) && ~isempty(sel_nonstm)
        selectivity_meta(j,16) = mean(sel_stim_leftALM(NdselRALM,cols1),'all','omitnan') / ...
                                 mean(sel_nonstm(NdselRALM,cols1),'all','omitnan');
    else
        selectivity_meta(j,16) = NaN;
    end

    % stimulus modularity
    if ~isempty(NdselLALM_stimulus) && ~isempty(sel_stimulus_right) && ~isempty(sel_stimulus_nostim)
        selectivity_meta(j,17) = mean(sel_stimulus_right(NdselLALM_stimulus,cols1),'all','omitnan') / ...
                                 mean(sel_stimulus_nostim(NdselLALM_stimulus,cols1),'all','omitnan');
    else
        selectivity_meta(j,17) = NaN;
    end
    if ~isempty(NdselRALM_stimulus) && ~isempty(sel_stimulus_left) && ~isempty(sel_stimulus_nostim)
        selectivity_meta(j,18) = mean(sel_stimulus_left(NdselRALM_stimulus,cols1),'all','omitnan') / ...
                                 mean(sel_stimulus_nostim(NdselRALM_stimulus,cols1),'all','omitnan');
    else
        selectivity_meta(j,18) = NaN;
    end

    % choice modularity
    if ~isempty(NdselLALM_choice) && ~isempty(sel_choice_right) && ~isempty(sel_choice_nostim)
        selectivity_meta(j,19) = mean(sel_choice_right(NdselLALM_choice,cols1),'all','omitnan') / ...
                                 mean(sel_choice_nostim(NdselLALM_choice,cols1),'all','omitnan');
    else
        selectivity_meta(j,19) = NaN;
    end
    if ~isempty(NdselRALM_choice) && ~isempty(sel_choice_left) && ~isempty(sel_choice_nostim)
        selectivity_meta(j,20) = mean(sel_choice_left(NdselRALM_choice,cols1),'all','omitnan') / ...
                                 mean(sel_choice_nostim(NdselRALM_choice,cols1),'all','omitnan');
    else
        selectivity_meta(j,20) = NaN;
    end

    % meta info
    mouse = obj2.Mice_session_all(j,1);
    session = obj2.Mice_session_all(j,2);
    session_str = num2str(session);
    obj3 = load([filepath,'\analyses_scripts\data\','meta_data_BAYLORGC',num2str(mouse),'_',session_str(1:4),'_',session_str(5:6),'_',session_str(7:8),'.mat']);
    selectivity_meta(j,23) = obj3.meta_data.craniotomy;
    mouse_session_tmp = obj2.Mice_session_all(find(obj2.Mice_session_all(:,1)==mouse),2);
    selectivity_meta(j,24) = obj3.meta_data.penetration(find(mouse_session_tmp==session));
    selectivity_meta(j,25) = obj3.meta_data.task;
end

save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat'], 'selectivity_meta_str', 'selectivity_meta');
%% save good delay selective sessions
clear;
filepath='C:\Users\XinHao\Desktop';
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);

Nselect_session_task1=find(selectivity_meta(:,11)>0.5&selectivity_meta(:,13)>0.5&selectivity_meta(:,12)>1&selectivity_meta(:,14)>1&...
    selectivity_meta(:,9)>=5&selectivity_meta(:,10)>=5&sum(selectivity_meta(:,9:10),2)>20&selectivity_meta(:,4)>0.6&selectivity_meta(:,5)>0.6&...
    selectivity_meta(:,24)==1&selectivity_meta(:,23)>0.5&selectivity_meta(:,25)==1);%49 sessions for standard task
Nselect_session_task2=find(selectivity_meta(:,11)>0.5&selectivity_meta(:,13)>0.5&selectivity_meta(:,12)>1&selectivity_meta(:,14)>1&...
    selectivity_meta(:,9)>=5&selectivity_meta(:,10)>=5&sum(selectivity_meta(:,9:10),2)>20&selectivity_meta(:,4)>0.6&selectivity_meta(:,5)>0.6&...
    selectivity_meta(:,24)==1&selectivity_meta(:,23)>0.5&selectivity_meta(:,25)==2);%6 sessions for reversed contingency task
Nselect_session_task3=find(selectivity_meta(:,11)>0.5&selectivity_meta(:,13)>0.5&selectivity_meta(:,12)>1&selectivity_meta(:,14)>1&...
    selectivity_meta(:,9)>=5&selectivity_meta(:,10)>=5&sum(selectivity_meta(:,9:10),2)>20&selectivity_meta(:,4)>0.6&selectivity_meta(:,5)>0.6&...
    selectivity_meta(:,24)==1&selectivity_meta(:,23)>0.5&selectivity_meta(:,25)==3);%30 sessions for fully reversed task
Nselect_session_task4=find(selectivity_meta(:,11)>0.5&selectivity_meta(:,13)>0.5&selectivity_meta(:,12)>1&selectivity_meta(:,14)>1&...
    selectivity_meta(:,9)>=5&selectivity_meta(:,10)>=5&sum(selectivity_meta(:,9:10),2)>20&selectivity_meta(:,4)>0.6&selectivity_meta(:,5)>0.6&...
    selectivity_meta(:,24)==1&selectivity_meta(:,23)>0.5&selectivity_meta(:,25)==4);%14 sessions for reversed tactile stimulus task
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);

%% 原文章的散点图
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);

figure;
for sub = 1:3
subplot(2,2,sub);
SiL=selectivity_meta(Nselect_session_task1, 2*sub+13);
SiR=selectivity_meta(Nselect_session_task1, 2*sub+14);
SiL(find(SiL>1))=1;%cap to between 0 and 1
SiR(find(SiR>1))=1;%cap to between 0 and 1
SiL(find(SiL<0))=0;%cap to between 0 and 1
SiR(find(SiR<0))=0;%cap to between 0 and 1
scatter(SiL,SiR,'ko');hold on; 
plot([0 1],[0 1],'k--');hold on;
 xlim([0 1]);
 ylim([0 1]);
xlabel('LALM modularity');
ylabel('RALM modularity');
end
saveas(gcf,'FigureScatter.fig','fig');


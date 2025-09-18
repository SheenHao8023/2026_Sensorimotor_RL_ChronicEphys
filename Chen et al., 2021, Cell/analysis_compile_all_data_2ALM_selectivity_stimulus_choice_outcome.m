%% establish, compute and save stimulus choice outcome selectivity variables
clear all
close all

addpath('.\helper_functions\')
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

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

    };

Stimulus_w=[];
Stimulus_dc_time=[];
Stimulus_dc_sig_time=[];
Stimulus_dc_sd_time=[];
Stimulus_sel_nonstm = [];

Choice_w=[];
Choice_dc_time=[];
Choice_dc_sig_time=[];
Choice_dc_sd_time=[];
Choice_sel_nonstm = [];

Outcome_w=[];
Outcome_dc_time=[];
Outcome_dc_sig_time=[];
Outcome_dc_sd_time=[];
Outcome_sel_nonstm = [];

sel_unitlist=[];
Mice_session_all=[];
Hemisphere_all=[];
sig_selective_all=[];
N_trials_all0=[];
FR_pref_all=[];

N_animals = length(Animals_list);  
n_session=0;

for I_animal = 1:N_animals
obj=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_PSTH_',Animals_list{I_animal},'.mat']);
Animals_list{I_animal}

for i_session=1:length(obj.PSTH_all)
    
    i_session
    n_session=n_session+1;
    
    N_trials_all=obj.N_trials_all{i_session};
    sig_selective=obj.Sig_selective_all{i_session};         %[sample delay resposne sdr]
    FR_pref=obj.FR_pref_all{i_session};
    Celltype=obj.CellType_all{i_session};
    PSTH=obj.PSTH_all{i_session};
    Hemisphere=obj.Hemisphere_all{i_session};
    Mouse_session=obj.Mice_all(i_session,:);
    t=obj.t;
    
for j=1:size(PSTH,1)
    for k=1:size(PSTH,2)
        if isnan(mean(mean(PSTH{j,k})))
            PSTH{j,k}=[];
        end
    end
end

sig_selective(isnan(sig_selective))=0;
i_selective = (sig_selective(:,1)|sig_selective(:,2)|sig_selective(:,3)) & Celltype'==1 & (abs(FR_pref(:,1))>0|abs(FR_pref(:,2))>0|abs(FR_pref(:,3))>0);
i_n_trial = sum(N_trials_all(1:size(PSTH,1),:),2)>100 & N_trials_all(1:size(PSTH,1),1)>5 & N_trials_all(1:size(PSTH,1),2)>5 & sum(N_trials_all(1:size(PSTH,1),[3 11]),2)>=2 & sum(N_trials_all(1:size(PSTH,1),[4 12]),2)>=2 & sum(N_trials_all(1:size(PSTH,1),[5 13]),2)>=2 & sum(N_trials_all(1:size(PSTH,1),[6 14]),2)>=2 & sum(N_trials_all(1:size(PSTH,1),[7 15]),2)>=2 & sum(N_trials_all(1:size(PSTH,1),[8 16]),2)>=2;%
n_useful=find(i_n_trial&i_selective)';



for i_cell =n_useful%
%     i_cell
    %both correct and error trials
    
    if ~isempty([PSTH{i_cell,1};PSTH{i_cell,9}])&~isempty([PSTH{i_cell,2};PSTH{i_cell,10}])%&...
               
        if size(PSTH{i_cell,9},1)>=2&size(PSTH{i_cell,10},1)>=2
            
            %for stimulus selectivity and cccp
            stimulus_w_time_tmp=[];%prefered direction, w>0 prefer right, w<0 prefer left
            stimulus_dc_time_tmp=[];%cccp combined conditions choice probability
            stimulus_dc_sig_time_tmp=[];%p value for cccp
            stimulus_dc_sd_time_tmp=[];%standard deviation for cccp
            stimulus_sel_nonstm_tmp = [];%selectivity
            
            for i_time=1:51%bin step 100ms, 51 bins
                tbin=100;
                psth1=[];
                psth2=[];
                psth3=[];
                psth4=[];
                
                if ~isempty(PSTH{i_cell,1})
                    psth1=nanmean(PSTH{i_cell,1}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,10})
                    psth2=nanmean(PSTH{i_cell,10}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,9})
                    psth3=nanmean(PSTH{i_cell,9}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,2})
                    psth4=nanmean(PSTH{i_cell,2}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                
                
                dpsth=[];
                for ktmp=1:10
                krand1=randsample(length(psth1),length(psth1),0);
                krand2=randsample(length(psth2),length(psth2),0);
                krand3=randsample(length(psth3),length(psth3),0);
                krand4=randsample(length(psth4),length(psth4),0);
                krand11=krand1(1:floor(length(krand1)/2));
                krand12=krand1(floor(length(krand1)/2)+1:end);
                krand21=krand2(1:floor(length(krand2)/2));
                krand22=krand2(floor(length(krand2)/2)+1:end);
                krand31=krand3(1:floor(length(krand3)/2));
                krand32=krand3(floor(length(krand3)/2)+1:end);
                krand41=krand4(1:floor(length(krand4)/2));
                krand42=krand4(floor(length(krand4)/2)+1:end);
                w1=sign(nanmean(psth1(krand11))-nanmean(psth2(krand21))+nanmean(psth3(krand31))-nanmean(psth4(krand41)));
                w2=sign(nanmean(psth1(krand12))-nanmean(psth2(krand22))+nanmean(psth3(krand32))-nanmean(psth4(krand42)));
                w=sign(nanmean(psth1)-nanmean(psth2)+nanmean(psth3)-nanmean(psth4));
                dpsth11=w2*(repmat(psth1(krand11),1,length(psth2(krand21)))-repmat(psth2(krand21)',length(psth1(krand11)),1));
                dpsth12=w2*(repmat(psth3(krand31),1,length(psth4(krand41)))-repmat(psth4(krand41)',length(psth3(krand31)),1));
                dpsth21=w1*(repmat(psth1(krand12),1,length(psth2(krand22)))-repmat(psth2(krand22)',length(psth1(krand12)),1));
                dpsth22=w1*(repmat(psth3(krand32),1,length(psth4(krand42)))-repmat(psth4(krand42)',length(psth3(krand32)),1));
                dpsth=[dpsth;reshape(dpsth11,prod(size(dpsth11)),1);reshape(dpsth12,prod(size(dpsth12)),1);...
                    reshape(dpsth21,prod(size(dpsth21)),1);reshape(dpsth22,prod(size(dpsth22)),1)];
                end
                dc_tmp=(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth);%cccp combined conditions choice probability of raw data
                
                
                stimulus_w_time_tmp(i_time)=w;
                for rep=1:100%change to 1000 will take very long time
                    psth_tmp=[psth1;psth2];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth1_tmp=psth_tmp(krand(1:length(psth1)));
                    psth2_tmp=psth_tmp(krand(length(psth1)+1:end));
                    psth_tmp=[psth3;psth4];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth3_tmp=psth_tmp(krand(1:length(psth3)));
                    psth4_tmp=psth_tmp(krand(length(psth3)+1:end));
                    dpsth1=(repmat(psth1_tmp,1,length(psth2_tmp))-repmat(psth2_tmp',length(psth1_tmp),1));
                    dpsth2=(repmat(psth3_tmp,1,length(psth4_tmp))-repmat(psth4_tmp',length(psth3_tmp),1));
                    dpsth=[reshape(dpsth1,prod(size(dpsth1)),1);reshape(dpsth2,prod(size(dpsth2)),1)];
                    dc_tmp=[dc_tmp;(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth)];%cccp combined conditions choice probability of shuffled data
                end
            
            stimulus_dc_time_tmp(i_time)=dc_tmp(1);%cccp combined conditions choice probability of raw data
            stimulus_dc_sig_time_tmp(i_time)=length(find(dc_tmp(2:end)>dc_tmp(1)))/rep;%p value
            stimulus_dc_sd_time_tmp(i_time)=nanstd(dc_tmp(2:end));
            end

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
                w1tmp=nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,10}(krand21,:),1)+nanmean(PSTH{i_cell,9}(krand31,:),1)-nanmean(PSTH{i_cell,2}(krand41,:),1);
                w1=sign(nanmean([nanmean(w1tmp(find(t>-3&t<-1.7))) nanmean(w1tmp(find(t>-1.7&t<0))) nanmean(w1tmp(find(t>0&t<1.5)))]));
                w2tmp=nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,10}(krand22,:),1)+nanmean(PSTH{i_cell,9}(krand32,:),1)-nanmean(PSTH{i_cell,2}(krand42,:),1);
                w2=sign(nanmean([nanmean(w2tmp(find(t>-3&t<-1.7))) nanmean(w2tmp(find(t>-1.7&t<0))) nanmean(w2tmp(find(t>0&t<1.5)))]));
                stimulus_sel_nonstm_tmp = [stimulus_sel_nonstm_tmp;w2*(nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,10}(krand21,:),1)+nanmean(PSTH{i_cell,9}(krand31,:),1)-nanmean(PSTH{i_cell,2}(krand41,:),1))/2;
                    w1*(nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,10}(krand22,:),1)+nanmean(PSTH{i_cell,9}(krand32,:),1)-nanmean(PSTH{i_cell,2}(krand42,:),1))/2];
            end
            stimulus_sel_nonstm_tmp=nanmean(stimulus_sel_nonstm_tmp,1);%stimulus selectivity
         
            %for choice selectivity and cccp
            choice_w_time_tmp=[];
            choice_dc_time_tmp=[];
            choice_dc_sig_time_tmp=[];
            choice_dc_sd_time_tmp=[];
            choice_sel_nonstm_tmp = [];
            for i_time=1:51%bin 100ms 51 bins
                tbin=100;
                psth1=[];
                psth2=[];
                psth3=[];
                psth4=[];
                if ~isempty(PSTH{i_cell,1})
                    psth1=nanmean(PSTH{i_cell,1}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,9})
                    psth2=nanmean(PSTH{i_cell,9}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,10})
                    psth3=nanmean(PSTH{i_cell,10}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,2})
                    psth4=nanmean(PSTH{i_cell,2}(:,[1:tbin]+(i_time-1)*tbin),2);
                end

                dpsth=[];
                for ktmp=1:10
                krand1=randsample(length(psth1),length(psth1),0);
                krand2=randsample(length(psth2),length(psth2),0);
                krand3=randsample(length(psth3),length(psth3),0);
                krand4=randsample(length(psth4),length(psth4),0);
                krand11=krand1(1:floor(length(krand1)/2));
                krand12=krand1(floor(length(krand1)/2)+1:end);
                krand21=krand2(1:floor(length(krand2)/2));
                krand22=krand2(floor(length(krand2)/2)+1:end);
                krand31=krand3(1:floor(length(krand3)/2));
                krand32=krand3(floor(length(krand3)/2)+1:end);
                krand41=krand4(1:floor(length(krand4)/2));
                krand42=krand4(floor(length(krand4)/2)+1:end);
                w1=sign(nanmean(psth1(krand11))-nanmean(psth2(krand21))+nanmean(psth3(krand31))-nanmean(psth4(krand41)));
                w2=sign(nanmean(psth1(krand12))-nanmean(psth2(krand22))+nanmean(psth3(krand32))-nanmean(psth4(krand42)));
                w=sign(nanmean(psth1)-nanmean(psth2)+nanmean(psth3)-nanmean(psth4));
                dpsth11=w2*(repmat(psth1(krand11),1,length(psth2(krand21)))-repmat(psth2(krand21)',length(psth1(krand11)),1));
                dpsth12=w2*(repmat(psth3(krand31),1,length(psth4(krand41)))-repmat(psth4(krand41)',length(psth3(krand31)),1));
                dpsth21=w1*(repmat(psth1(krand12),1,length(psth2(krand22)))-repmat(psth2(krand22)',length(psth1(krand12)),1));
                dpsth22=w1*(repmat(psth3(krand32),1,length(psth4(krand42)))-repmat(psth4(krand42)',length(psth3(krand32)),1));
                dpsth=[dpsth;reshape(dpsth11,prod(size(dpsth11)),1);reshape(dpsth12,prod(size(dpsth12)),1);...
                    reshape(dpsth21,prod(size(dpsth21)),1);reshape(dpsth22,prod(size(dpsth22)),1)];
                end
                dc_tmp=(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth);
                
                choice_w_time_tmp(i_time)=w;
                for rep=1:100%change to 1000 will take very long time
                    psth_tmp=[psth1;psth2];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth1_tmp=psth_tmp(krand(1:length(psth1)));
                    psth2_tmp=psth_tmp(krand(length(psth1)+1:end));
                    psth_tmp=[psth3;psth4];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth3_tmp=psth_tmp(krand(1:length(psth3)));
                    psth4_tmp=psth_tmp(krand(length(psth3)+1:end));
                    dpsth1=(repmat(psth1_tmp,1,length(psth2_tmp))-repmat(psth2_tmp',length(psth1_tmp),1));
                    dpsth2=(repmat(psth3_tmp,1,length(psth4_tmp))-repmat(psth4_tmp',length(psth3_tmp),1));
                    dpsth=[reshape(dpsth1,prod(size(dpsth1)),1);reshape(dpsth2,prod(size(dpsth2)),1)];
                    dc_tmp=[dc_tmp;(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth)];
                end
            choice_dc_time_tmp(i_time)=dc_tmp(1);
            choice_dc_sig_time_tmp(i_time)=length(find(dc_tmp(2:end)>dc_tmp(1)))/rep;
            choice_dc_sd_time_tmp(i_time)=nanstd(dc_tmp(2:end));
            end

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
                w1tmp=nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,9}(krand21,:),1)+nanmean(PSTH{i_cell,10}(krand31,:),1)-nanmean(PSTH{i_cell,2}(krand41,:),1);
                w1=sign(nanmean([nanmean(w1tmp(find(t>-3&t<-1.7))) nanmean(w1tmp(find(t>-1.7&t<0))) nanmean(w1tmp(find(t>0&t<1.5)))]));
                w2tmp=nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,9}(krand22,:),1)+nanmean(PSTH{i_cell,10}(krand32,:),1)-nanmean(PSTH{i_cell,2}(krand42,:),1);
                w2=sign(nanmean([nanmean(w2tmp(find(t>-3&t<-1.7))) nanmean(w2tmp(find(t>-1.7&t<0))) nanmean(w2tmp(find(t>0&t<1.5)))]));
                choice_sel_nonstm_tmp = [choice_sel_nonstm_tmp;w2*(nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,9}(krand21,:),1)+nanmean(PSTH{i_cell,10}(krand31,:),1)-nanmean(PSTH{i_cell,2}(krand41,:),1))/2;
                    w1*(nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,9}(krand22,:),1)+nanmean(PSTH{i_cell,10}(krand32,:),1)-nanmean(PSTH{i_cell,2}(krand42,:),1))/2];
            end
            
            choice_sel_nonstm_tmp =nanmean(choice_sel_nonstm_tmp,1);

            %for outcome selectivity and cccp
            outcome_w_time_tmp=[];
            outcome_dc_time_tmp=[];
            outcome_dc_sig_time_tmp=[];
            outcome_dc_sd_time_tmp=[];
            outcome_sel_nonstm_tmp = [];
            for i_time=1:51%bin 100ms 51 bins
                tbin=100;
                psth1=[];
                psth2=[];
                psth3=[];
                psth4=[];
                if ~isempty(PSTH{i_cell,1})
                    psth1=nanmean(PSTH{i_cell,1}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,9})
                    psth2=nanmean(PSTH{i_cell,9}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,2})
                    psth3=nanmean(PSTH{i_cell,2}(:,[1:tbin]+(i_time-1)*tbin),2);
                end
                if ~isempty(PSTH{i_cell,10})
                    psth4=nanmean(PSTH{i_cell,10}(:,[1:tbin]+(i_time-1)*tbin),2);
                end

                dpsth=[];
                for ktmp=1:10
                krand1=randsample(length(psth1),length(psth1),0);
                krand2=randsample(length(psth2),length(psth2),0);
                krand3=randsample(length(psth3),length(psth3),0);
                krand4=randsample(length(psth4),length(psth4),0);
                krand11=krand1(1:floor(length(krand2)/2));
                krand12=krand1(end-floor(length(krand2)/2):end);
                krand21=krand2(1:floor(length(krand2)/2));
                krand22=krand2(floor(length(krand2)/2)+1:end);
                krand31=krand3(1:floor(length(krand4)/2));
                krand32=krand3(end-floor(length(krand4)/2):end);
                krand41=krand4(1:floor(length(krand4)/2));
                krand42=krand4(floor(length(krand4)/2)+1:end);
                w1=sign(nanmean(psth1(krand11))-nanmean(psth2(krand21))+nanmean(psth3(krand31))-nanmean(psth4(krand41)));
                w2=sign(nanmean(psth1(krand12))-nanmean(psth2(krand22))+nanmean(psth3(krand32))-nanmean(psth4(krand42)));
                w=sign(nanmean(psth1)-nanmean(psth2)+nanmean(psth3)-nanmean(psth4));
                dpsth11=w2*(repmat(psth1(krand11),1,length(psth2(krand21)))-repmat(psth2(krand21)',length(psth1(krand11)),1));
                dpsth12=w2*(repmat(psth3(krand31),1,length(psth4(krand41)))-repmat(psth4(krand41)',length(psth3(krand31)),1));
                dpsth21=w1*(repmat(psth1(krand12),1,length(psth2(krand22)))-repmat(psth2(krand22)',length(psth1(krand12)),1));
                dpsth22=w1*(repmat(psth3(krand32),1,length(psth4(krand42)))-repmat(psth4(krand42)',length(psth3(krand32)),1));
                dpsth=[dpsth;reshape(dpsth11,prod(size(dpsth11)),1);reshape(dpsth12,prod(size(dpsth12)),1);...
                    reshape(dpsth21,prod(size(dpsth21)),1);reshape(dpsth22,prod(size(dpsth22)),1)];
                end
                dc_tmp=(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth);
                
                outcome_w_time_tmp(i_time)=w;
                for rep=1:100%change to 1000 will take very long time
                    psth_tmp=[psth1;psth2];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth1_tmp=psth_tmp(krand(1:length(psth1)));
                    psth2_tmp=psth_tmp(krand(length(psth1)+1:end));
                    psth_tmp=[psth3;psth4];
                    krand=randsample(length(psth_tmp),length(psth_tmp),0);
                    psth3_tmp=psth_tmp(krand(1:length(psth3)));
                    psth4_tmp=psth_tmp(krand(length(psth3)+1:end));
                    dpsth1=(repmat(psth1_tmp,1,length(psth2_tmp))-repmat(psth2_tmp',length(psth1_tmp),1));
                    dpsth2=(repmat(psth3_tmp,1,length(psth4_tmp))-repmat(psth4_tmp',length(psth3_tmp),1));
                    dpsth=[reshape(dpsth1,prod(size(dpsth1)),1);reshape(dpsth2,prod(size(dpsth2)),1)];
                    dc_tmp=[dc_tmp;(length(find(dpsth>0))+length(find(dpsth==0))/2)/length(dpsth)];
                end
            outcome_dc_time_tmp(i_time)=dc_tmp(1);
            outcome_dc_sig_time_tmp(i_time)=length(find(dc_tmp(2:end)>dc_tmp(1)))/rep;
            outcome_dc_sd_time_tmp(i_time)=nanstd(dc_tmp(2:end));
            end

            for ktmp=1:10
                krand=randsample(size(PSTH{i_cell,1},1),size(PSTH{i_cell,1},1),0);
                krand11=krand(1:floor(length(krand)/2));
                krand12=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,9},1),size(PSTH{i_cell,9},1),0);
                krand21=krand(1:floor(length(krand)/2));
                krand22=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,2},1),size(PSTH{i_cell,2},1),0);
                krand31=krand(1:floor(length(krand)/2));
                krand32=krand(floor(length(krand)/2)+1:end);
                krand=randsample(size(PSTH{i_cell,10},1),size(PSTH{i_cell,10},1),0);
                krand41=krand(1:floor(length(krand)/2));
                krand42=krand(floor(length(krand)/2)+1:end);
                w1tmp=nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,9}(krand21,:),1)+nanmean(PSTH{i_cell,2}(krand31,:),1)-nanmean(PSTH{i_cell,10}(krand41,:),1);
                w1=sign(nanmean([nanmean(w1tmp(find(t>-3&t<-1.7))) nanmean(w1tmp(find(t>-1.7&t<0))) nanmean(w1tmp(find(t>0&t<1.5)))]));
                w2tmp=nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,9}(krand22,:),1)+nanmean(PSTH{i_cell,2}(krand32,:),1)-nanmean(PSTH{i_cell,10}(krand42,:),1);
                w2=sign(nanmean([nanmean(w2tmp(find(t>-3&t<-1.7))) nanmean(w2tmp(find(t>-1.7&t<0))) nanmean(w2tmp(find(t>0&t<1.5)))]));
                outcome_sel_nonstm_tmp = [outcome_sel_nonstm_tmp;w2*(nanmean(PSTH{i_cell,1}(krand11,:),1)-nanmean(PSTH{i_cell,9}(krand21,:),1)+nanmean(PSTH{i_cell,2}(krand31,:),1)-nanmean(PSTH{i_cell,10}(krand41,:),1))/2;
                    w1*(nanmean(PSTH{i_cell,1}(krand12,:),1)-nanmean(PSTH{i_cell,9}(krand22,:),1)+nanmean(PSTH{i_cell,2}(krand32,:),1)-nanmean(PSTH{i_cell,10}(krand42,:),1))/2];
            end
            
            outcome_sel_nonstm_tmp =nanmean(outcome_sel_nonstm_tmp,1);
                
            Stimulus_dc_time(end+1,:)=stimulus_dc_time_tmp;
            Stimulus_dc_sig_time(end+1,:)=stimulus_dc_sig_time_tmp;
            Stimulus_dc_sd_time(end+1,:)=stimulus_dc_sd_time_tmp;
            Stimulus_sel_nonstm(end+1,:) = stimulus_sel_nonstm_tmp;
            Choice_dc_time(end+1,:)=choice_dc_time_tmp;
            Choice_dc_sig_time(end+1,:)=choice_dc_sig_time_tmp;
            Choice_dc_sd_time(end+1,:)=choice_dc_sd_time_tmp;
            Choice_sel_nonstm(end+1,:) =choice_sel_nonstm_tmp;
            Outcome_dc_time(end+1,:)=outcome_dc_time_tmp;
            Outcome_dc_sig_time(end+1,:)=outcome_dc_sig_time_tmp;
            Outcome_dc_sd_time(end+1,:)=outcome_dc_sd_time_tmp;
            Outcome_sel_nonstm(end+1,:) = outcome_sel_nonstm_tmp;
            sel_unitlist(end+1,:)=i_cell;
            Stimulus_w(end+1,:)=stimulus_w_time_tmp;
            Choice_w(end+1,:)=choice_w_time_tmp;
            Outcome_w(end+1,:)=outcome_w_time_tmp;
            Mice_session_all(end+1,:)=Mouse_session;
            Hemisphere_all(end+1,:)=Hemisphere(i_cell);
            sig_selective_all(end+1,:)=sig_selective(i_cell,:);
            N_trials_all0(end+1,:)=N_trials_all(i_cell,:);
            FR_pref_all(end+1,:)=FR_pref(i_cell,:);
        end
end
end
end
clear obj;
end
N_trials_all=N_trials_all0;
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_stimulus_choice_outcome.mat'],...
    't','sig_selective_all','N_trials_all','FR_pref_all','sel_unitlist','Stimulus_w','Choice_w','Outcome_w',...
    'Stimulus_dc_time','Stimulus_dc_sig_time','Stimulus_dc_sd_time','Stimulus_sel_nonstm',...
    'Choice_dc_time','Choice_dc_sig_time','Choice_dc_sd_time','Choice_sel_nonstm',...
    'Outcome_dc_time','Outcome_dc_sig_time','Outcome_dc_sd_time','Outcome_sel_nonstm',...
    'Mice_session_all','Hemisphere_all','-v7.3');

%% Figure 1C bottom, Figure S1D
clear;
filepath='C:\Users\XinHao\Desktop';

load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_stimulus_choice_outcome.mat']);
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);
%find high selective sessions as saved in selectivity_meta.mat
sessionall=ismember(Mice_session_all(:,1),selectivity_meta(Nselect_session_task1,1))&ismember(Mice_session_all(:,2),selectivity_meta(Nselect_session_task1,2));
figure;
for z=1:2
idx=find(sessionall==1&Hemisphere_all==z);
Stimulus_dc_time_tmp=Stimulus_dc_time(idx,:);
Stimulus_dc_sig_time_tmp=Stimulus_dc_sig_time(idx,:);
Stimulus_dc_sd_time_tmp=Stimulus_dc_sd_time(idx,:);
Stimulus_sel_nonstm_tmp = Stimulus_sel_nonstm(idx,:);
Choice_dc_time_tmp=Choice_dc_time(idx,:);
Choice_dc_sig_time_tmp=Choice_dc_sig_time(idx,:);
Choice_dc_sd_time_tmp=Choice_dc_sd_time(idx,:);
Choice_sel_nonstm_tmp=Choice_sel_nonstm(idx,:);
Outcome_dc_time_tmp=Outcome_dc_time(idx,:);
Outcome_dc_sig_time_tmp=Outcome_dc_sig_time(idx,:);
Outcome_dc_sd_time_tmp=Outcome_dc_sd_time(idx,:);
Outcome_sel_nonstm_tmp= Outcome_sel_nonstm(idx,:);
Stimulus_w_time_tmp=Stimulus_w(idx,:);
Choice_w_time_tmp=Choice_w(idx,:);
Outcome_w_time_tmp=Outcome_w(idx,:);

%
pvalue=0.05;
xt=([1:size(Stimulus_dc_time_tmp,2)]-abs(t(1))/0.1-1)*0.1;
tsd=[8:34];%sample delay, sd
tr=[35:length(xt)];%response, r
Stimulus_neuron_sd=find(mean(Stimulus_dc_sig_time_tmp(:,tsd),2)<pvalue&mean(Choice_dc_sig_time_tmp(:,tsd),2)>=pvalue&mean(Outcome_dc_sig_time_tmp(:,tr),2)>=pvalue);%stimulus selective only
Choice_neuron_sd=find(mean(Stimulus_dc_sig_time_tmp(:,tsd),2)>=pvalue&mean(Choice_dc_sig_time_tmp(:,tsd),2)<pvalue&mean(Outcome_dc_sig_time_tmp(:,tr),2)>=pvalue);%choice selective only
Outcome_neuron_r=find(mean(Stimulus_dc_sig_time_tmp(:,tsd),2)>=pvalue&mean(Choice_dc_sig_time_tmp(:,tsd),2)>=pvalue&mean(Outcome_dc_sig_time_tmp(:,tr),2)<pvalue);%outcome selective only

Stimulus_neuron=[];
Choice_neuron=[];
Outcome_neuron=[];
for k=1:length(xt)
Stimulus_neuron(:,k)=(Stimulus_dc_sig_time_tmp(:,k)<pvalue&Choice_dc_sig_time_tmp(:,k)>=pvalue&Outcome_dc_sig_time_tmp(:,k)>=pvalue);
Choice_neuron(:,k)=(Stimulus_dc_sig_time_tmp(:,k)>=pvalue&Choice_dc_sig_time_tmp(:,k)<pvalue&Outcome_dc_sig_time_tmp(:,k)>=pvalue);
Outcome_neuron(:,k)=(Stimulus_dc_sig_time_tmp(:,k)>=pvalue&Choice_dc_sig_time_tmp(:,k)>=pvalue&Outcome_dc_sig_time_tmp(:,k)<pvalue);
end

for j=1:3
    switch j
        case 1
            neuron_sdr=Stimulus_neuron_sd;
            neuron=Stimulus_neuron;
            dc_sig_time=Stimulus_dc_sig_time_tmp;
            dc_time=Stimulus_dc_time_tmp;
            w_time=Stimulus_w_time_tmp;
            sel_nonstm=Stimulus_sel_nonstm_tmp;
            tit='stimulus coding';
        case 2
            neuron_sdr=Choice_neuron_sd;
            neuron=Choice_neuron;
            dc_sig_time=Choice_dc_sig_time_tmp;
            dc_time=Choice_dc_time_tmp;
            w_time=Choice_w_time_tmp;
            sel_nonstm=Choice_sel_nonstm_tmp;
            tit='choice coding';
        case 3
            neuron_sdr=Outcome_neuron_r;
            neuron=Outcome_neuron;
            dc_sig_time=Outcome_dc_sig_time_tmp;
            dc_time=Outcome_dc_time_tmp;
            w_time=Outcome_w_time_tmp;
            sel_nonstm=Outcome_sel_nonstm_tmp;
            tit='outcome coding';
    end


subplot(2,3,z);hold on;
sig=[];
for k=1:length(xt)
sig(k)=length(find(neuron(:,k)))/length(neuron);
end
if j==1
plot(xt,sig,'g');hold on;%stimulus selective
end
if j==2
plot(xt,sig,'b');hold on;%choice selective
end
line([0 0],[-0.3 0.3],'color','k')
line([-1.7 -1.7],[-0.3 0.3],'color','k')
line([-3 -3],[-0.3 0.3],'color','k')
plot([-3.5 2],[0.05 0.05],'k--')
xlim([-3.5 2]);
ylim([0 0.3]);
xlabel('time(s)');
ylabel('frac. sig. neuron');
if z==1
    title('Left ALM');
else
    title('Right ALM');
end

subplot(2,3,3+j);hold on;
if z==1
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'k');
outlinebounds(hl,hp);
hold on;
else
    [hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'b');
outlinebounds(hl,hp);
hold on;
end
Ylim=round(1.3*max([max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlim([-3.5 1.5]);
xlabel('time(s)');
ylabel('selectivity');
if j==1
title('stimulus selectivity');
else if j==2
        title('choice selectivity');
    else
        title('outcome selectivity');
    end
end

end
end
saveas(gcf,'Figure1Cbottom_Figure S1D.fig','fig');
%% FigureS1E
clear;
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_stimulus_choice_outcome.mat']);
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);
%compute stimulus/choice selectivity for each session
LALM_Stimulus_sel_nonstm=[];
LALM_Choice_sel_nonstm=[];

RALM_Stimulus_sel_nonstm=[];
RALM_Choice_sel_nonstm=[];

for p=Nselect_session_task1'

idx=find(ismember(Mice_session_all(:,1),selectivity_meta(p,1))&ismember(Mice_session_all(:,2),selectivity_meta(p,2))&Hemisphere_all==1);
LALM_Stimulus_sel_nonstm=[LALM_Stimulus_sel_nonstm;mean(Stimulus_sel_nonstm(idx,:))];
LALM_Choice_sel_nonstm=[LALM_Choice_sel_nonstm;mean(Choice_sel_nonstm(idx,:))];
idx=find(ismember(Mice_session_all(:,1),selectivity_meta(p,1))&ismember(Mice_session_all(:,2),selectivity_meta(p,2))&Hemisphere_all==2);
RALM_Stimulus_sel_nonstm=[RALM_Stimulus_sel_nonstm;mean(Stimulus_sel_nonstm(idx,:))];
RALM_Choice_sel_nonstm=[RALM_Choice_sel_nonstm;mean(Choice_sel_nonstm(idx,:))];

end
%
figure;
for j=1:3
    switch j
        case 1
            twin=find(t>-3&t<-1.7);
            tit='sample';
        case 2
            twin=find(t>-1.7&t<0);
            tit='delay';
        case 3
            twin=find(t>0&t<1.5);
            tit='response';
    end
subplot(2,3,j);
scatter(nanmean(LALM_Stimulus_sel_nonstm(:,twin),2),nanmean(RALM_Stimulus_sel_nonstm(:,twin),2));hold on;
plot([0 2],[0 2],'k--');hold on;
[h,p]=ttest(nanmean(LALM_Stimulus_sel_nonstm(:,twin),2),nanmean(RALM_Stimulus_sel_nonstm(:,twin),2));
text(0.5,0.5,['p:',num2str(p)]);
xlabel('Left ALM StiSel');
ylabel('Right ALM StiSel');
title(tit);
if j==1
    xlim([0 1]);
    ylim([0 1]);
else if j==2
    xlim([0 2]);
    ylim([0 2]);
    else
    xlim([0 2]);
    ylim([0 2]);
    end
end
subplot(2,3,3+j);
scatter(nanmean(LALM_Choice_sel_nonstm(:,twin),2),nanmean(RALM_Choice_sel_nonstm(:,twin),2));hold on;
plot([0 2],[0 2],'k--');hold on;
[h,p]=ttest(nanmean(LALM_Choice_sel_nonstm(:,twin),2),nanmean(RALM_Choice_sel_nonstm(:,twin),2));
text(0.5,0.5,['p:',num2str(p)]);
xlabel('Left ALM ChoiceSel');
ylabel('Right ALM ChoiceSel');
title(tit);
if j==1
    xlim([0 1]);
    ylim([0 1]);
else if j==2
    xlim([0 2]);
    ylim([0 2]);
    else
    xlim([0 3]);
    ylim([0 3]);
    end
end
end
saveas(gcf,'FigureS1E.fig','fig');


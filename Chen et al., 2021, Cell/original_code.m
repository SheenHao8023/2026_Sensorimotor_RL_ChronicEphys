% establish, compute and save trial type selectivity variables
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
%
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

NdselLALM_all={};
NdselRALM_all={};
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
%opto tagging, identify putative inhibitory neurons strongly excited by light
tsti=find(t>-1.7&t<-0.9);%early delay
optotag_cell=zeros(size(PSTH,1),1);
for j=1:size(PSTH,1)
psthnostim=[PSTH{j,1};PSTH{j,2};PSTH{j,9};PSTH{j,10}];
psthleft=[PSTH{j,3};PSTH{j,4};PSTH{j,11};PSTH{j,12}];
psthright=[PSTH{j,5};PSTH{j,6};PSTH{j,13};PSTH{j,14}];
psthbilat=[PSTH{j,7};PSTH{j,8};PSTH{j,15};PSTH{j,16}];

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

N=size(PSTH,1);

%compute psth for control, stimleft, stimright trials
Psth_nostim=[];
Psth_stim_left=[];
Psth_stim_right=[];
Psth_nostim_mean=[];
Psth_stim_left_mean=[];
Psth_stim_right_mean=[];
Psth_useful_unit_hemisphere=[];

k=0;
for j=n_useful      
        psthnostim=[PSTH{j,1};PSTH{j,2};PSTH{j,9};PSTH{j,10}];
        psthleft=[PSTH{j,3};PSTH{j,4};PSTH{j,11};PSTH{j,12}];
        psthright=[PSTH{j,5};PSTH{j,6};PSTH{j,13};PSTH{j,14}];  
        psthbilat=[PSTH{j,7};PSTH{j,8};PSTH{j,15};PSTH{j,16}];
        if ~isempty(psthnostim)&~isempty(psthleft)&~isempty(psthright)&~isempty(psthbilat)
        k=k+1;
        Psth_nostim(k,:)=nanmean(psthnostim);
        Psth_stim_left(k,:)=nanmean(psthleft);
        Psth_stim_right(k,:)=nanmean(psthright);
        
        Psth_nostim_mean(k)=nanmean(nanmean(psthnostim(:,tsti)));
        Psth_stim_left_mean(k)=nanmean(nanmean(psthleft(:,tsti)));
        Psth_stim_right_mean(k)=nanmean(nanmean(psthright(:,tsti)));
        
        Psth_useful_unit_hemisphere(k)=Hemisphere(j);
        end
end

Psth_nostim_all{n_session}=Psth_nostim;
Psth_stim_left_all{n_session}=Psth_stim_left;
Psth_stim_right_all{n_session}=Psth_stim_right;
Psth_nostim_mean_all{n_session}=Psth_nostim_mean;
Psth_stim_left_mean_all{n_session}=Psth_stim_left_mean;
Psth_stim_right_mean_all{n_session}=Psth_stim_right_mean;
Psth_useful_unit_hemisphere_all{n_session}=Psth_useful_unit_hemisphere;


%compute selectivity
sel_nonstm = [];
sel_stim_leftALM = [];
sel_stim_rightALM = [];
sel_stim_Bilat = [];
sig_selective_Selected = [];
sel_unitlist=[];
sel_unit_hemisphere=[];

for i_cell =intersect(n_useful,find(i_selective & i_n_trial)')%use loose selective excitatory neurons
    
    if ~isempty([PSTH{i_cell,1};PSTH{i_cell,9}])&~isempty([PSTH{i_cell,2};PSTH{i_cell,10}])&...
            ~isempty([PSTH{i_cell,3};PSTH{i_cell,11}])&~isempty([PSTH{i_cell,4};PSTH{i_cell,12}])&...
            ~isempty([PSTH{i_cell,5};PSTH{i_cell,13}])&~isempty([PSTH{i_cell,6};PSTH{i_cell,14}])%&...
            
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

end
clear obj
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_',Animals_list{I_animal},'.mat'],...
    'Psth_nostim_all','Psth_stim_left_all','Psth_stim_right_all','Psth_nostim_mean_all',...
    'Psth_stim_left_mean_all','Psth_stim_right_mean_all','Psth_useful_unit_hemisphere_all',...
    'sel_nonstm_all','sel_stim_leftALM_all','sel_stim_rightALM_all','sel_stim_Bilat_all',...
    'sig_selective_Selected_all','sel_unitlist_all','sel_unit_hemisphere_all',...
    'NdselLALM_all','NdselRALM_all','NsdrselLALM_all','NsdrselRALM_all',...
    'Mice_session_all','t','tsti');
clearvars -except Animals_list I_animal N_animals filepath
end
%% combine all sessions selectivity related variables and save them to one file
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

NdselLALM_all={};
NdselRALM_all={};
NsdrselLALM_all={};
NsdrselRALM_all={};
Mice_session_all=[];

n_session=0;
N_animals = length(Animals_list); 

for I_animal = 1:N_animals
Animals_list{I_animal}
obj=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_',Animals_list{I_animal},'.mat']);

for i_session=1:length(obj.Psth_nostim_all)
    i_session
    n_session=n_session+1;
    Psth_nostim_all{n_session}=obj.Psth_nostim_all{i_session};
    Psth_stim_left_all{n_session}=obj.Psth_stim_left_all{i_session};
    Psth_stim_right_all{n_session}=obj.Psth_stim_right_all{i_session};
    Psth_nostim_mean_all{n_session}=obj.Psth_nostim_mean_all{i_session};
    Psth_stim_left_mean_all{n_session}=obj.Psth_stim_left_mean_all{i_session};
    Psth_stim_right_mean_all{n_session}=obj.Psth_stim_right_mean_all{i_session};
    Psth_useful_unit_hemisphere_all{n_session}=obj.Psth_useful_unit_hemisphere_all{i_session};
    
    sel_nonstm_all{n_session}=obj.sel_nonstm_all{i_session};
    sel_stim_leftALM_all{n_session}=obj.sel_stim_leftALM_all{i_session};
    sel_stim_rightALM_all{n_session}=obj.sel_stim_rightALM_all{i_session};
    sel_stim_Bilat_all{n_session}=obj.sel_stim_Bilat_all{i_session};
    sig_selective_Selected_all{n_session}=obj.sig_selective_Selected_all{i_session};
    sel_unitlist_all{n_session}=obj.sel_unitlist_all{i_session};
    sel_unit_hemisphere_all{n_session}=obj.sel_unit_hemisphere_all{i_session};
    
    NdselLALM_all{n_session}=obj.NdselLALM_all{i_session};
    NdselRALM_all{n_session}=obj.NdselRALM_all{i_session};
    NsdrselLALM_all{n_session}=obj.NsdrselLALM_all{i_session};
    NsdrselRALM_all{n_session}=obj.NsdrselRALM_all{i_session};
    Mice_session_all(n_session,:)=obj.Mice_session_all(i_session,:);
end
end
t=obj.t;
tsti=obj.tsti;
clear obj
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat'],...
    'Psth_nostim_all','Psth_stim_left_all','Psth_stim_right_all','Psth_nostim_mean_all',...
    'Psth_stim_left_mean_all','Psth_stim_right_mean_all','Psth_useful_unit_hemisphere_all',...
    'sel_nonstm_all','sel_stim_leftALM_all','sel_stim_rightALM_all','sel_stim_Bilat_all',...
    'sig_selective_Selected_all','sel_unitlist_all','sel_unit_hemisphere_all',...
    'NdselLALM_all','NdselRALM_all','NsdrselLALM_all','NsdrselRALM_all',...
    'Mice_session_all','t','tsti');
%% delay selectivity, modularity, robustness, and other meta information
clear;
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

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
twin3=[-0.1 0];
obj1=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_behavior.mat']);
obj2=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
t=obj2.t;
selectivity_meta=[];
for j=1:size(obj2.Mice_session_all,1)
    j
    selectivity_meta(j,1:2)=obj2.Mice_session_all(j,1:2);%mice session name
    if obj1.behavior_performance_allSession(j,1)==obj2.Mice_session_all(j,1)
        if obj1.behavior_performance_allSession(j,2)==obj2.Mice_session_all(j,2)
           selectivity_meta(j,3:8)=obj1.behavior_performance_allSession(j,3:end);% behavior performance
        end%'perf_ctr','perf_ctr_Rtrial','perf_ctr_Ltrial','perf_stiLeft','perf_stiRight','perf_stiBilat'
    end
    NdselLALM=obj2.NdselLALM_all{j};
    NdselRALM=obj2.NdselRALM_all{j};
    sel_nonstm=obj2.sel_nonstm_all{j};
    sel_stim_leftALM=obj2.sel_stim_leftALM_all{j};
    sel_stim_rightALM=obj2.sel_stim_rightALM_all{j};
    selectivity_meta(j,9:10)=[length(NdselLALM) length(NdselRALM)];%number of delay selective neurons 'num_dselLALM','num_dselRALM'
    selectivity_meta(j,11)=mean(mean(sel_nonstm(NdselLALM,find(t>twin1(1)&t<twin1(2)))));%LALM early delay 'LALM_earlyD_sel'
    selectivity_meta(j,12)=mean(mean(sel_nonstm(NdselLALM,find(t>twin2(1)&t<twin2(2)))));%LALM late delay 'LALM_lateD_sel'
    selectivity_meta(j,13)=mean(mean(sel_nonstm(NdselRALM,find(t>twin1(1)&t<twin1(2)))));%RALM early delay 'RALM_earlyD_sel'
    selectivity_meta(j,14)=mean(mean(sel_nonstm(NdselRALM,find(t>twin2(1)&t<twin2(2)))));%RALM late delay 'RALM_lateD_sel'
    
    selectivity_meta(j,15)=mean(mean(sel_stim_rightALM(NdselLALM,find(t>twin1(1)&t<twin1(2)))))/mean(mean(sel_nonstm(NdselLALM,find(t>twin1(1)&t<twin1(2)))));%'LALM_modularity'
    selectivity_meta(j,16)=mean(mean(sel_stim_leftALM(NdselRALM,find(t>twin1(1)&t<twin1(2)))))/mean(mean(sel_nonstm(NdselRALM,find(t>twin1(1)&t<twin1(2)))));%'RALM_modularity'
    
    selectivity_meta(j,17)=mean(mean(sel_stim_leftALM(NdselLALM,find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm(NdselLALM,find(t>twin3(1)&t<twin3(2)))));%'LALMstiLSelrec'
    selectivity_meta(j,18)=mean(mean(sel_stim_rightALM(NdselLALM,find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm(NdselLALM,find(t>twin3(1)&t<twin3(2)))));%'LALMstiRSelrec'

    selectivity_meta(j,19)=mean(mean(sel_stim_leftALM(NdselRALM,find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm(NdselRALM,find(t>twin3(1)&t<twin3(2)))));%'RALMstiLSelrec'
    selectivity_meta(j,20)=mean(mean(sel_stim_rightALM(NdselRALM,find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm(NdselRALM,find(t>twin3(1)&t<twin3(2)))));%'RALMstiRSelrec'

    selectivity_meta(j,21)=mean(mean(sel_stim_leftALM([NdselLALM;NdselRALM],find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm([NdselLALM;NdselRALM],find(t>twin3(1)&t<twin3(2)))));%'LRALMstiLSelrec'
    selectivity_meta(j,22)=mean(mean(sel_stim_rightALM([NdselLALM;NdselRALM],find(t>twin3(1)&t<twin3(2)))))/mean(mean(sel_nonstm([NdselLALM;NdselRALM],find(t>twin3(1)&t<twin3(2)))));%'LRALMstiRSelrec'

    mouse=obj2.Mice_session_all(j,1);
    session=obj2.Mice_session_all(j,2);
    session_str=num2str(session);
    obj3=load([filepath,'\data\','meta_data_BAYLORGC',num2str(mouse),'_',session_str(1:4),'_',session_str(5:6),'_',session_str(7:8),'.mat']);
    selectivity_meta(j,23)=obj3.meta_data.craniotomy;%'craniotomy_quality'
    mouse_session_tmp=obj2.Mice_session_all(find(obj2.Mice_session_all(:,1)==mouse),2);
    selectivity_meta(j,24)=obj3.meta_data.penetration(find(mouse_session_tmp==session));%'penetration_quality'
    selectivity_meta(j,25)=obj3.meta_data.task;%'task_type'
    selectivity_meta(j,26)=length(obj3.meta_data.animalStrain);%'strain' 1 VGAT, 2 PVReachR
    if strcmp(obj3.meta_data.sex,'f')
    selectivity_meta(j,27)=1;%'sex' 1 female, 2 male
    else
    selectivity_meta(j,27)=2;%'sex' 1 female, 2 male
    end
    selectivity_meta(j,28)=obj3.meta_data.age;%'age'
    selectivity_meta(j,29)=obj3.meta_data.trainingDays;%'trainingdays'
    
    tsti=find(t>-1.6&t<-1);%cut 0.1 s at the two ends of [-1.7:-0.9]s considering the 0.2s average window
    tmp=mean(obj2.Psth_nostim_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==1&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2);%Left ALM control early delay rate
    selectivity_meta(j,30)=nanmean(mean(obj2.Psth_stim_left_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==1&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2)./tmp);%Left ALM sti left early delay rate
    selectivity_meta(j,31)=nanmean(mean(obj2.Psth_stim_right_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==1&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2)./tmp);%Left ALM sti right early delay rate
    tmp=mean(obj2.Psth_nostim_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==2&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2);%Right ALM control early delay rate
    selectivity_meta(j,32)=nanmean(mean(obj2.Psth_stim_left_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==2&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2)./tmp);%Right ALM sti left early delay rate
    selectivity_meta(j,33)=nanmean(mean(obj2.Psth_stim_right_all{j}(find(obj2.Psth_useful_unit_hemisphere_all{j}==2&mean(obj2.Psth_nostim_all{j}(:,find(t<-3)),2)'>=1),tsti),2)./tmp);%Right ALM sti right early delay rate
end
save([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat'],'selectivity_meta_str','selectivity_meta');
%% save good delay selective sessions
clear;
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

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
%% plot Figure 2B C  
clear;
filepath='C:\Users\XinHao\Desktop';

load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
figure;
for m=1:2
    if m==1
       N=47 %modular session 47 BAYLORGC86_2019_07_31; asym session 35 BAYLORGC50_2018_12_07
    else
        N=35
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
subplot(2,3,(m-1)*3+j);
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
else if j==2
ylabel('RALM selectivity(spikes/s)');
    end   
end
Nsel=size(sel_nonstm,1)
if j==1
    title(['StiRight',' Nsel=',num2str(Nsel)]);
end
if j==2
    title(['StiLeft',' Nsel=',num2str(Nsel)]);
end
Si=[];
for k=1:1000
    krand=randsample(size(sel_sti,1),size(sel_sti,1),1);
    Si(k)=mean(mean(sel_sti(krand,find(t>-1.7&t<-0.9))))/mean(mean(sel_nonstm(krand,find(t>-1.7&t<-0.9))));
end
subplot(2,3,(m-1)*3+3);hold on;
bar(j,mean(Si));hold on;
errorbar(j,mean(Si),std(Si));hold on;
ylabel('modularity');
xlabel('Left Right (ALM)');
end
end
saveas(gcf,'Figure2B_C.fig','fig');
%% plot Figure 3, scatter plots could be slightly different to the figure in the paper due to random sampling during selectivity calculation
% each dot position may be slightly different for every run of selectivity
% calculation due to random split sampling of trials for each neuron (see
% the 1st part code)
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);

figure;
for j=1:4
    switch j
        case 1
            N=Nselect_session_task1;
        case 2
            N=Nselect_session_task3;
        case 3
            N=Nselect_session_task2;
        case 4
            N=Nselect_session_task4;
    end  
subplot(2,2,j);
SiL=selectivity_meta(N,15);
SiR=selectivity_meta(N,16);
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
saveas(gcf,'Figure3.fig','fig');
%% plot figure 4 B E, they could be slightly different to the figure in the paper due to random sampling during selectivity calculation
% standard task
Group=Nselect_session_task1;
SiL=selectivity_meta(:,15);%left ALM modularity
SiR=selectivity_meta(:,16);%right ALM modularity
SiL(find(SiL>1))=1;%cap to between 0 and 1
SiR(find(SiR>1))=1;%cap to between 0 and 1
SiL(find(SiL<0))=0;%cap to between 0 and 1
SiR(find(SiR<0))=0;%cap to between 0 and 1
group1=Group(find(SiL(Group)>0.7));
group2=Group(find(SiR(Group)>0.47));
asymmetry=SiL(group1)-SiR(group1);
modularity=SiL(group2)+SiR(group2);
[Sa,Ia]=sort(asymmetry);
[Sm,Im]=sort(modularity);
group1=group1(Ia(end-6:end));%top 7 asymmetric sessions
group2=group2(Im(end-6:end));%top 7 modular sessions

figure;
for G=1:2
    if G==1
        group=group2;
    else if G==2
            group=group1;
        end
    end
    LALM_sel_nonstm=[];
    LALM_sel_stiL=[];
    RALM_sel_nonstm=[];
    RALM_sel_stiL=[];
    for g=group'
        LALM_sel_nonstm=[LALM_sel_nonstm;sel_nonstm_all{g}(NdselLALM_all{g},:)];
        LALM_sel_stiL=[LALM_sel_stiL;sel_stim_leftALM_all{g}(NdselLALM_all{g},:)];
        RALM_sel_nonstm=[RALM_sel_nonstm;sel_nonstm_all{g}(NdselRALM_all{g},:)];
        RALM_sel_stiL=[RALM_sel_stiL;sel_stim_leftALM_all{g}(NdselRALM_all{g},:)];
    end

subplot(2,3,(G-1)*3+1);
[hl,hp]=boundedline(t,mean(LALM_sel_nonstm),std(LALM_sel_nonstm)/sqrt(size(LALM_sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(LALM_sel_stiL),std(LALM_sel_stiL)/sqrt(size(LALM_sel_stiL,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(LALM_sel_stiL)) max(mean(LALM_sel_nonstm))]));
ylim([-1, 4]);
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('LALM selectivity(spikes/s)');

subplot(2,3,(G-1)*3+2);
[hl,hp]=boundedline(t,mean(RALM_sel_nonstm),std(RALM_sel_nonstm)/sqrt(size(RALM_sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(RALM_sel_stiL),std(RALM_sel_stiL)/sqrt(size(RALM_sel_stiL,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(RALM_sel_stiL)) max(mean(RALM_sel_nonstm))]));
ylim([-1, 4]);
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('RALM selectivity(spikes/s)');

subplot(2,3,(G-1)*3+3);
bar(1,mean(SiR(group)));hold on;
scatter(ones(1,length(group)),SiR(group),'ko');hold on;
ylim([0 1]);
xlabel('Right ALM');
ylabel('modularity');
end
saveas(gcf,'Figure4BE.fig','fig');
%% Figure 4G-I, Figure S3H-I
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);
xdata1=selectivity_meta(Nselect_session_task1,21);%LRALM stiL selectivity recovery
xdata2=selectivity_meta(Nselect_session_task1,22);%LRALM stiR selectivity recovery
xdata3=selectivity_meta(Nselect_session_task3,21);
xdata4=selectivity_meta(Nselect_session_task3,22);
ydata1=selectivity_meta(Nselect_session_task1,6)./selectivity_meta(Nselect_session_task1,3);%performance stiLeft recovery
ydata2=selectivity_meta(Nselect_session_task1,7)./selectivity_meta(Nselect_session_task1,3);%performance stiRight recovery
ydata3=selectivity_meta(Nselect_session_task3,6)./selectivity_meta(Nselect_session_task3,3);
ydata4=selectivity_meta(Nselect_session_task3,7)./selectivity_meta(Nselect_session_task3,3);
zdata1=selectivity_meta(Nselect_session_task1,16);%right ALM modularity
zdata2=selectivity_meta(Nselect_session_task1,15);%left ALM modularity
zdata3=selectivity_meta(Nselect_session_task3,16);
zdata4=selectivity_meta(Nselect_session_task3,15);
xdata1(find(xdata1<0))=0;
xdata1(find(xdata1>1))=1;
xdata2(find(xdata2<0))=0;
xdata2(find(xdata2>1))=1;
xdata3(find(xdata3<0))=0;
xdata3(find(xdata3>1))=1;
xdata4(find(xdata4<0))=0;
xdata4(find(xdata4>1))=1;
ydata1(find(ydata1<0))=0;
ydata1(find(ydata1>1))=1;
ydata2(find(ydata2<0))=0;
ydata2(find(ydata2>1))=1;
ydata3(find(ydata3<0))=0;
ydata3(find(ydata3>1))=1;
ydata4(find(ydata4<0))=0;
ydata4(find(ydata4>1))=1;
zdata1(find(zdata1<0))=0;
zdata1(find(zdata1>1))=1;
zdata2(find(zdata2<0))=0;
zdata2(find(zdata2>1))=1;
zdata3(find(zdata3<0))=0;
zdata3(find(zdata3>1))=1;
zdata4(find(zdata4<0))=0;
zdata4(find(zdata4>1))=1;
xdata=[xdata1;xdata2;xdata3;xdata4];
ydata=[ydata1;ydata2;ydata3;ydata4];
zdata=[zdata1;zdata2;zdata3;zdata4];

%normalzied rate during photoinhibition
rdata1=selectivity_meta(Nselect_session_task1,30);%LALM stiL
rdata2=selectivity_meta(Nselect_session_task1,31);%LALM stiR
rdata3=selectivity_meta(Nselect_session_task3,30);%LALM stiL
rdata4=selectivity_meta(Nselect_session_task3,31);%LALM stiR
rdata5=selectivity_meta(Nselect_session_task1,32);%RALM stiL
rdata6=selectivity_meta(Nselect_session_task1,33);%RALM stiR
rdata7=selectivity_meta(Nselect_session_task3,32);%RALM stiL
rdata8=selectivity_meta(Nselect_session_task3,33);%RALM stiR

%Figure 4G-I
figure;
subplot(1,3,1);
scatter(xdata1,ydata1,'ko','filled');hold on;
scatter(xdata2,ydata2,'ko');hold on;
scatter(xdata3,ydata3,'k^','filled');hold on;
scatter(xdata4,ydata4,'k^');hold on;
FM = fitlm(xdata',ydata','linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
text(0,0.8,['R2:',num2str(FM.Rsquared.Ordinary)]);hold on;
[r,p]=corr(xdata,ydata,'type','Pearson');
text(0,0.7,['r:',num2str(r)]);hold on;
text(0,0.6,['p:',num2str(p)]);hold on;
xlim([0 1]);
ylim([0.5 1]);
xlabel('neuronal recovery');
ylabel('behavior recovery');
subplot(1,3,2);
scatter(zdata1,xdata1,'ko','filled');hold on;
scatter(zdata2,xdata2,'ko');hold on;
scatter(zdata3,xdata3,'k^','filled');hold on;
scatter(zdata4,xdata4,'k^');hold on;
FM = fitlm(zdata',xdata','linear');
plot(zdata,FM.Coefficients{2,'Estimate'}*zdata+FM.Coefficients{1,'Estimate'},'k');hold on;
text(0,0.7,['R2:',num2str(FM.Rsquared.Ordinary)]);hold on;
[r,p]=corr(zdata,xdata,'type','Pearson');
text(0,0.5,['r:',num2str(r)]);hold on;
text(0,0.3,['p:',num2str(p)]);hold on;
xlim([0 1]);
ylim([0 1]);
xlabel('modularity');
ylabel('neuronal recovery');
subplot(1,3,3);
scatter(zdata1,ydata1,'ko','filled');hold on;
scatter(zdata2,ydata2,'ko');hold on;
scatter(zdata3,ydata3,'k^','filled');hold on;
scatter(zdata4,ydata4,'k^');hold on;
FM = fitlm(zdata',ydata','linear');
plot(zdata,FM.Coefficients{2,'Estimate'}*zdata+FM.Coefficients{1,'Estimate'},'k');hold on;
text(0,0.8,['R2:',num2str(FM.Rsquared.Ordinary)]);hold on;
[r,p]=corr(zdata,ydata,'type','Pearson');
text(0,0.7,['r:',num2str(r)]);hold on;
text(0,0.6,['p:',num2str(p)]);hold on;
xlim([0 1]);
ylim([0.5 1]);
xlabel('modularity');
ylabel('behavior recovery');
saveas(gcf,'Figure4G_I.fig','fig');

%Figure S3H-I
figure;
subplot(1,4,1);
scatter([zdata1;zdata3],[rdata1;rdata3],'ko');hold on;
plot([0 1],[1 1],'k--');hold on;
xlim([0 1]);
ylim([0 1.5]);
xlabel('RALM modularity');
ylabel('Norm. spike rate');
title('LALM stim left');
subplot(1,4,2);
scatter([zdata1;zdata3],[rdata5;rdata7],'ko');hold on;
plot([0 1],[1 1],'k--');hold on;
xlim([0 1]);
ylim([0 1.5]);
xlabel('RALM modularity');
ylabel('Norm. spike rate');
title('RALM stim left');
subplot(1,4,3);
scatter([zdata2;zdata4],[rdata2;rdata4],'ko');hold on;
plot([0 1],[1 1],'k--');hold on;
xlim([0 1]);
ylim([0 1.5]);
xlabel('LALM modularity');
ylabel('Norm. spike rate');
title('LALM stim right');
subplot(1,4,4);
scatter([zdata2;zdata4],[rdata6;rdata8],'ko');hold on;
plot([0 1],[1 1],'k--');hold on;
xlim([0 1]);
ylim([0 1.5]);
xlabel('LALM modularity');
ylabel('Norm. spike rate');
title('RALM stim right');
saveas(gcf,'FigureS3H_I.fig','fig');

%% FigureS2G FigureS3B-C FigureS3D-E
clear all;
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

obj=load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity.mat']);
load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);
t=obj.t;

%FigureS2G
group=[Nselect_session_task1;Nselect_session_task3]';
LALM_sel_nonstm=[];
LALM_sel_stim_Bilat=[];
RALM_sel_nonstm=[];
RALM_sel_stim_Bilat=[];
for j=group
    LALM_sel_nonstm=[LALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselLALM_all{j},:)];
    LALM_sel_stim_Bilat=[LALM_sel_stim_Bilat;obj.sel_stim_Bilat_all{j}(obj.NdselLALM_all{j},:)];
    RALM_sel_nonstm=[RALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselRALM_all{j},:)];
    RALM_sel_stim_Bilat=[RALM_sel_stim_Bilat;obj.sel_stim_Bilat_all{j}(obj.NdselRALM_all{j},:)];      
end
for k=1:2
    if k==1
        sel_nonstm=LALM_sel_nonstm;
        sel_stim_Bilat=LALM_sel_stim_Bilat;
        tit='LALM stim Bilat';
    else
        sel_nonstm=RALM_sel_nonstm;
        sel_stim_Bilat=RALM_sel_stim_Bilat;
        tit='RALM stim Bilat';
    end
subplot(1,2,k);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_stim_Bilat),std(sel_stim_Bilat)/sqrt(size(sel_stim_Bilat,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_stim_Bilat)) max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('selec(spks/s)');
xlim([-3.5 2])
ylim([-1 Ylim])
line([-1.7 -0.9],[0.8 0.8]*Ylim,'color',[0 0.6 1],'linewidth',3)
title(tit);
end
saveas(gcf,'FigureS2G.fig','fig');

SiL=selectivity_meta(:,15);%left ALM modularity
SiR=selectivity_meta(:,16);%right ALM modularity
SiL(find(SiL>1))=1;%cap to between 0 and 1
SiR(find(SiR>1))=1;%cap to between 0 and 1
SiL(find(SiL<0))=0;%cap to between 0 and 1
SiR(find(SiR<0))=0;%cap to between 0 and 1
group=[Nselect_session_task1;Nselect_session_task3]';

%FigureS3B-C
figure;
[S,I]=sort(SiR(group));
group2=group(I(1:floor(length(I)/4)-1));
group1=group(I(round(length(I)/4*3)+1:end));
group2x=group(I(floor(length(I)/4)));
group1x=group(I(round(length(I)/4*3)));

subplot(3,3,1:3);
Si=[0:0.05:1];
Sc=hist(SiR(group),Si);
bar(Si,Sc);hold on;
plot([SiR(group2x) SiR(group2x)],[0 5]);hold on;
plot([SiR(group1x) SiR(group1x)],[0 5]);hold on;
xlabel('RALM modularity');
ylabel('sessions');

subplot(3,3,4);
bar(1,mean(SiR(group1)));hold on;
scatter(ones(1,length(group1)),SiR(group1),'ko');hold on;
ylim([0 1]);
xlabel('RALM');
ylabel('modularity');
title(['n=',num2str(length(group1))]);

LALM_sel_nonstm=[];
LALM_sel_stim_left=[];
RALM_sel_nonstm=[];
RALM_sel_stim_left=[];
for j=group1
    LALM_sel_nonstm=[LALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselLALM_all{j},:)];
    LALM_sel_stim_left=[LALM_sel_stim_left;obj.sel_stim_leftALM_all{j}(obj.NdselLALM_all{j},:)];
    RALM_sel_nonstm=[RALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselRALM_all{j},:)];
    RALM_sel_stim_left=[RALM_sel_stim_left;obj.sel_stim_leftALM_all{j}(obj.NdselRALM_all{j},:)];      
end
for k=1:2
    if k==1
        sel_nonstm=LALM_sel_nonstm;
        sel_stim_left=LALM_sel_stim_left;
        tit='LALM stim left';
    else
        sel_nonstm=RALM_sel_nonstm;
        sel_stim_left=RALM_sel_stim_left;
        tit='RALM stim left';
    end
subplot(3,3,4+k);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_stim_left),std(sel_stim_left)/sqrt(size(sel_stim_left,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_stim_left)) max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('selec(spks/s)');
xlim([-3.5 2])
ylim([-1 Ylim])
line([-1.7 -0.9],[0.8 0.8]*Ylim,'color',[0 0.6 1],'linewidth',3)
title(tit);
end

subplot(3,3,7);
bar(1,mean(SiR(group2)));hold on;
scatter(ones(1,length(group2)),SiR(group2),'ko');hold on;
ylim([0 1]);
xlabel('RALM');
ylabel('modularity');
title(['n=',num2str(length(group2))]);

LALM_sel_nonstm=[];
LALM_sel_stim_left=[];
RALM_sel_nonstm=[];
RALM_sel_stim_left=[];
for j=group2
    LALM_sel_nonstm=[LALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselLALM_all{j},:)];
    LALM_sel_stim_left=[LALM_sel_stim_left;obj.sel_stim_leftALM_all{j}(obj.NdselLALM_all{j},:)];
    RALM_sel_nonstm=[RALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselRALM_all{j},:)];
    RALM_sel_stim_left=[RALM_sel_stim_left;obj.sel_stim_leftALM_all{j}(obj.NdselRALM_all{j},:)];      
end
for k=1:2
    if k==1
        sel_nonstm=LALM_sel_nonstm;
        sel_stim_left=LALM_sel_stim_left;
        tit='LALM stim left';
    else
        sel_nonstm=RALM_sel_nonstm;
        sel_stim_left=RALM_sel_stim_left;
        tit='RALM stim left';
    end
subplot(3,3,7+k);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_stim_left),std(sel_stim_left)/sqrt(size(sel_stim_left,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_stim_left)) max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('selec(spks/s)');
xlim([-3.5 2])
ylim([-1 Ylim])
line([-1.7 -0.9],[0.8 0.8]*Ylim,'color',[0 0.6 1],'linewidth',3)
title(tit);
end
saveas(gcf,'FigureS3BC.fig','fig');

%FigureS3D-E
figure;
[S,I]=sort(SiL(group));
group2=group(I(1:floor(length(I)/4)-1));
group1=group(I(round(length(I)/4*3)+1:end));
group2x=group(I(floor(length(I)/4)));
group1x=group(I(round(length(I)/4*3)));

subplot(3,3,1:3);
Si=[0:0.05:1];
Sc=hist(SiL(group),Si);
bar(Si,Sc);hold on;
plot([SiL(group2x) SiL(group2x)],[0 5]);hold on;
plot([SiL(group1x) SiL(group1x)],[0 5]);hold on;
xlabel('LALM modularity');
ylabel('sessions');

subplot(3,3,4);
bar(1,mean(SiL(group1)));hold on;
scatter(ones(1,length(group1)),SiL(group1),'ko');hold on;
ylim([0 1]);
xlabel('LALM');
ylabel('modularity');
title(['n=',num2str(length(group1))]);

LALM_sel_nonstm=[];
LALM_sel_stim_right=[];
RALM_sel_nonstm=[];
RALM_sel_stim_right=[];
for j=group1
    LALM_sel_nonstm=[LALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselLALM_all{j},:)];
    LALM_sel_stim_right=[LALM_sel_stim_right;obj.sel_stim_rightALM_all{j}(obj.NdselLALM_all{j},:)];
    RALM_sel_nonstm=[RALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselRALM_all{j},:)];
    RALM_sel_stim_right=[RALM_sel_stim_right;obj.sel_stim_rightALM_all{j}(obj.NdselRALM_all{j},:)];      
end
for k=1:2
    if k==1
        sel_nonstm=LALM_sel_nonstm;
        sel_stim_right=LALM_sel_stim_right;
        tit='LALM stim right';
    else
        sel_nonstm=RALM_sel_nonstm;
        sel_stim_right=RALM_sel_stim_right;
        tit='RALM stim right';
    end
subplot(3,3,4+k);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_stim_right),std(sel_stim_right)/sqrt(size(sel_stim_right,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_stim_right)) max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('selec(spks/s)');
xlim([-3.5 2])
ylim([-1 Ylim])
line([-1.7 -0.9],[0.8 0.8]*Ylim,'color',[0 0.6 1],'linewidth',3)
title(tit);
end

subplot(3,3,7);
bar(1,mean(SiL(group2)));hold on;
scatter(ones(1,length(group2)),SiL(group2),'ko');hold on;
ylim([0 1]);
xlabel('LALM');
ylabel('modularity');
title(['n=',num2str(length(group2))]);

LALM_sel_nonstm=[];
LALM_sel_stim_right=[];
RALM_sel_nonstm=[];
RALM_sel_stim_right=[];
for j=group2
    LALM_sel_nonstm=[LALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselLALM_all{j},:)];
    LALM_sel_stim_right=[LALM_sel_stim_right;obj.sel_stim_rightALM_all{j}(obj.NdselLALM_all{j},:)];
    RALM_sel_nonstm=[RALM_sel_nonstm;obj.sel_nonstm_all{j}(obj.NdselRALM_all{j},:)];
    RALM_sel_stim_right=[RALM_sel_stim_right;obj.sel_stim_rightALM_all{j}(obj.NdselRALM_all{j},:)];      
end
for k=1:2
    if k==1
        sel_nonstm=LALM_sel_nonstm;
        sel_stim_right=LALM_sel_stim_right;
        tit='LALM stim right';
    else
        sel_nonstm=RALM_sel_nonstm;
        sel_stim_right=RALM_sel_stim_right;
        tit='RALM stim right';
    end
subplot(3,3,7+k);
[hl,hp]=boundedline(t,mean(sel_nonstm),std(sel_nonstm)/sqrt(size(sel_nonstm,1)),'-k','transparency',0.1);
outlinebounds(hl,hp);
hold on;
[hl,hp]=boundedline(t,mean(sel_stim_right),std(sel_stim_right)/sqrt(size(sel_stim_right,1)),'-b','transparency',0.1);
outlinebounds(hl,hp);
hold on;
Ylim=round(1.3*max([max(mean(sel_stim_right)) max(mean(sel_nonstm))]));
line([0 0],[-1 Ylim],'color','k')
line([-1.7 -1.7],[-1 Ylim],'color','k')
line([-3 -3],[-1 Ylim],'color','k')
line([min(t) max(t)],[0 0],'color','k')
xlabel('time (s)')
ylabel('selec(spks/s)');
xlim([-3.5 2])
ylim([-1 Ylim])
line([-1.7 -0.9],[0.8 0.8]*Ylim,'color',[0 0.6 1],'linewidth',3)
title(tit);
end
saveas(gcf,'FigureS3DE.fig','fig');

%% FigureS5A-E  FigureS5F-G
clear;
filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021';

load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_selectivity_meta.mat']);
%filepath='Y:\Lab data\Guang Chen\ChenKangEtAl_2021'; set this in case the
%old filepath in *_meta.mat is different to the current one

load([filepath,'\analyses_scripts\analysis_compile_all_data_2ALM_CDprojection_CC.mat']);
N=length(Nselect_session_task1);

mouseID=selectivity_meta(Nselect_session_task1,1);
SiL=selectivity_meta(Nselect_session_task1,15);%left ALM modularity
SiR=selectivity_meta(Nselect_session_task1,16);%right ALM modularity
SiL(find(SiL>1))=1;%cap to between 0 and 1
SiR(find(SiR>1))=1;%cap to between 0 and 1
SiL(find(SiL<0))=0;%cap to between 0 and 1
SiR(find(SiR<0))=0;%cap to between 0 and 1
strain=selectivity_meta(Nselect_session_task1,26);
gender=selectivity_meta(Nselect_session_task1,27);
age=selectivity_meta(Nselect_session_task1,28);
trainday=selectivity_meta(Nselect_session_task1,29);

%FigureS5A-E
LRALMstiLSelrec=selectivity_meta(Nselect_session_task1,21);
LRALMstiRSelrec=selectivity_meta(Nselect_session_task1,22);
stiLPerfrec=selectivity_meta(Nselect_session_task1,6)./selectivity_meta(Nselect_session_task1,3);
stiRPerfrec=selectivity_meta(Nselect_session_task1,7)./selectivity_meta(Nselect_session_task1,3);
LRALMstiLSelrec(find(LRALMstiLSelrec<0))=0;
LRALMstiLSelrec(find(LRALMstiLSelrec>1))=1;
LRALMstiRSelrec(find(LRALMstiRSelrec<0))=0;
LRALMstiRSelrec(find(LRALMstiRSelrec>1))=1;
stiLPerfrec(find(stiLPerfrec<0))=0;
stiLPerfrec(find(stiLPerfrec>1))=1;
stiRPerfrec(find(stiRPerfrec<0))=0;
stiRPerfrec(find(stiRPerfrec>1))=1;

rankCC=squeeze(LRALM_CDproj_cc_t(1:N,1,3,end));%standard task, control, both trial type, late delay CDprojection, left right ALM correlation

SiL_crossmice=[];
SiR_crossmice=[];
SiL_withinmice=[];
SiR_withinmice=[];

LRALMstiLSelrec_crossmice=[];
LRALMstiRSelrec_crossmice=[];
LRALMstiLSelrec_withinmice=[];
LRALMstiRSelrec_withinmice=[];

stiLPerfrec_crossmice=[];
stiRPerfrec_crossmice=[];
stiLPerfrec_withinmice=[];
stiRPerfrec_withinmice=[];

rankCC_crossmice=[];
trainday_crossmice=[];

mice=unique(mouseID);
for j=1:length(mice)
    idx=find(mouseID==mice(j));
    SiL_crossmice(j)=mean(SiL(idx));
    SiR_crossmice(j)=mean(SiR(idx));
    SiL_withinmice=[SiL_withinmice;SiL(idx)-mean(SiL(idx))];
    SiR_withinmice=[SiR_withinmice;SiR(idx)-mean(SiR(idx))];

    LRALMstiLSelrec_crossmice(j)=mean(LRALMstiLSelrec(idx));
    LRALMstiRSelrec_crossmice(j)=mean(LRALMstiRSelrec(idx));
    LRALMstiLSelrec_withinmice=[LRALMstiLSelrec_withinmice;LRALMstiLSelrec(idx)-mean(LRALMstiLSelrec(idx))];
    LRALMstiRSelrec_withinmice=[LRALMstiRSelrec_withinmice;LRALMstiRSelrec(idx)-mean(LRALMstiRSelrec(idx))];
    
    stiLPerfrec_crossmice(j)=mean(stiLPerfrec(idx));
    stiRPerfrec_crossmice(j)=mean(stiRPerfrec(idx));
    stiLPerfrec_withinmice=[stiLPerfrec_withinmice;stiLPerfrec(idx)-mean(stiLPerfrec(idx))];
    stiRPerfrec_withinmice=[stiRPerfrec_withinmice;stiRPerfrec(idx)-mean(stiRPerfrec(idx))];
      
    rankCC_crossmice(j)=mean(rankCC(idx));
    trainday_crossmice(j)=mean(trainday(idx));
end

figure;
subplot(3,3,2);
scatter(SiL,SiR,'k');hold on;
plot([0 1],[0 1],'k--');hold on;
xlim([0 1]);
ylim([0 1]);
xlabel('modularity, LeftALM');
ylabel('modularity, RightALM');
subplot(3,3,3);
scatter(SiL_withinmice,SiR_withinmice,'k');hold on;
plot([-0.5 0.5],[-0.5 0.5],'k--');hold on;
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
xlabel('delta modularity across sessions, LeftALM');
ylabel('delta modularity across sessions, RightALM');
subplot(3,3,4);
xdata=[SiL_crossmice SiR_crossmice];
ydata=[LRALMstiRSelrec_crossmice LRALMstiLSelrec_crossmice];
scatter(SiL_crossmice,LRALMstiRSelrec_crossmice,'ko');hold on;
scatter(SiR_crossmice,LRALMstiLSelrec_crossmice,'ko','filled');hold on;
FM = fitlm(xdata,ydata,'linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata',ydata','type','Pearson');
text(0.6,0.2,['r:',num2str(r)]);hold on;
text(0.6,0.1,['p:',num2str(p)]);hold on;
xlim([0 1]);
ylim([0 1]);
xlabel('modularity');
ylabel('neuronal recovery');
subplot(3,3,5);
xdata=[SiL_crossmice SiR_crossmice];
ydata=[stiRPerfrec_crossmice stiLPerfrec_crossmice];
scatter(SiL_crossmice,stiRPerfrec_crossmice,'ko');hold on;
scatter(SiR_crossmice,stiLPerfrec_crossmice,'ko','filled');hold on;
FM = fitlm(xdata,ydata,'linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata',ydata','type','Pearson');
text(0.6,0.2,['r:',num2str(r)]);hold on;
text(0.6,0.1,['p:',num2str(p)]);hold on;
xlim([0 1]);
ylim([0 1]);
xlabel('modularity');
ylabel('behavior recovery');
subplot(3,3,6);
xdata=SiL_crossmice-SiR_crossmice;
ydata=rankCC_crossmice;
scatter(xdata,ydata,'ko');hold on;
FM = fitlm(xdata,ydata,'linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata',ydata','type','Pearson');
text(0.6,0.2,['r:',num2str(r)]);hold on;
text(0.6,0.1,['p:',num2str(p)]);hold on;
xlim([-0.2 0.6]);
ylim([0 1]);
xlabel('asymmetry');
ylabel('rank correlation');

idx=find(stiRPerfrec_withinmice~=0);%only use mice with >=2 sessions
subplot(3,3,7);
xdata=[SiL_withinmice(idx);SiR_withinmice(idx)];
ydata=[LRALMstiRSelrec_withinmice(idx);LRALMstiLSelrec_withinmice(idx)];
scatter(SiL_withinmice(idx),LRALMstiRSelrec_withinmice(idx),'ko');hold on;
scatter(SiR_withinmice(idx),LRALMstiLSelrec_withinmice(idx),'ko','filled');hold on;
FM = fitlm(xdata',ydata','linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata,ydata,'type','Pearson');
text(0.2,-0.2,['r:',num2str(r)]);hold on;
text(0.2,-0.4,['p:',num2str(p)]);hold on;
xlim([-0.4 0.4]);
ylim([-0.4 0.4]);
xlabel('delta modularity');
ylabel('delta neuronal recovery');
subplot(3,3,8);
xdata=[SiL_withinmice(idx);SiR_withinmice(idx)];
ydata=[stiRPerfrec_withinmice(idx);stiLPerfrec_withinmice(idx)];
scatter(SiL_withinmice(idx),stiRPerfrec_withinmice(idx),'ko');hold on;
scatter(SiR_withinmice(idx),stiLPerfrec_withinmice(idx),'ko','filled');hold on;
FM = fitlm(xdata',ydata','linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata,ydata,'type','Pearson');
text(0.3,-0.1,['r:',num2str(r)]);hold on;
text(0.3,-0.2,['p:',num2str(p)]);hold on;
xlim([-0.4 0.4]);
ylim([-0.2 0.2]);
xlabel('delta modularity');
ylabel('delta behavior recovery');


asymmetry=SiL-SiR;
asymmetrysort=sort(asymmetry);
asym1=asymmetrysort(floor(length(asymmetrysort)/3));
asym2=asymmetrysort(round(length(asymmetrysort)/3*2)+1);

asymmetry_withinmice=[];
rankCC_withinmice=[];
for j=1:length(asymmetry)
    if length(find(mouseID==mouseID(j)))>1&max(asymmetry(find(mouseID==mouseID(j))))>asym2%mouse with at least two sessions and at least one high asymmetric session
       asymmetry_withinmice(j)=(asymmetry(j)-mean(asymmetry(find(mouseID==mouseID(j)))));
       rankCC_withinmice(j)=(rankCC(j)-mean(rankCC(find(mouseID==mouseID(j)))));      
    else
       asymmetry_withinmice(j)=NaN;
       rankCC_withinmice(j)=NaN;
    end
end

subplot(3,3,9);
xdata=asymmetry_withinmice(find(~isnan(asymmetry_withinmice)));
ydata=rankCC_withinmice(find(~isnan(asymmetry_withinmice)));
scatter(xdata,ydata,'ko');hold on;
FM = fitlm(xdata,ydata,'linear');
plot(xdata,FM.Coefficients{2,'Estimate'}*xdata+FM.Coefficients{1,'Estimate'},'k');hold on;
[r,p]=corr(xdata',ydata','type','Pearson');
text(0.2,-0.1,['r:',num2str(r)]);hold on;
text(0.2,-0.2,['p:',num2str(p)]);hold on;
xlim([-0.4 0.4]);
ylim([-0.3 0.3]);
xlabel('delta asymmetry');
ylabel('delta rank correlation');
saveas(gcf,'FigureS5A_E.fig','fig');


%FigureS5F-G
figure;
for j=1:2
    if j==1
        ydata=SiL;
        fname='Modularity,leftALM';
    else
        ydata=SiR;
        fname='Modularity,rightALM';
    end
subplot(2,4,(j-1)*4+1);
bar([1 2],[mean(ydata(find(strain==1))) mean(ydata(find(strain==2)))]);hold on;
scatter(strain,ydata);
[h,p]=ttest2(ydata(find(strain==1)),ydata(find(strain==2)));%
text(1,0.5,['p=',num2str(p)]);
ylabel(fname);
xlabel('VGAT/PVReachR');
subplot(2,4,(j-1)*4+2);
bar([1 2],[mean(ydata(find(gender==1))) mean(ydata(find(gender==2)))]);hold on;
scatter(gender,ydata);
[h,p]=ttest2(ydata(find(gender==1)),ydata(find(gender==2)));%
text(1,0.5,['p=',num2str(p)]);
ylabel(fname);
xlabel('female/male');
subplot(2,4,(j-1)*4+3);
scatter(age,ydata,'ko');hold on;
xlim([100 400]);
ylim([0 1]);
ylabel(fname);
xlabel('Age(days)');
end

for j=1:2
    if j==1
        ydata=SiL_crossmice;
        fname='Modularity,leftALM';
    else
        ydata=SiR_crossmice;
        fname='Modularity,rightALM';
    end
subplot(2,4,(j-1)*4+4);
scatter(trainday_crossmice,ydata,'ko');hold on;
xlim([0 80]);
ylim([0 1]);
ylabel(fname);
xlabel('training days to perf70%');
end
saveas(gcf,'FigureS5F_G.fig','fig');

function func_process_voltage_traces_kilosort(input_filelist, output_dir, ch_spike, trigger_ch, bitcode_ch, probe, shank, laser_trial)

    rng(1);   
    tag = round(rand(1)*1000); % 生成的文件名包含随机数，前面固定了随机种子，所以不会跑一次生成不同的文件

    % check input
    if size(input_filelist,1)~=1 & size(input_filelist,2)~=1
        error('input_filelist must be a vector')
    end
    if size(input_filelist,2)~=1
        input_filelist = input_filelist';
    end

    % sort the files by number
    input_fileNo_tmp = [];
    for i_file_tmp = 1:size(input_filelist,1)
        filename_tmp = input_filelist(i_file_tmp).name;
        i_str1 = strfind(filename_tmp,'_');
        i_str2 = strfind(filename_tmp,'.rhd');
        input_fileNo_tmp(i_file_tmp) = str2num([filename_tmp(i_str1(end-1)+1:i_str1(end)-1) filename_tmp(i_str1(end)+1:i_str2(1)-1)]);
    end
    [file_num_sorted_tmp, i_file_sorted_tmp] = sort(input_fileNo_tmp);
    if sum(diff(file_num_sorted_tmp)>1)>0
        warning('files are not sequencial, file numbers will not be guaranteed when bitcode fails')
    end
    input_filelist = input_filelist(i_file_sorted_tmp);    
    clear *_tmp

    % start to process data
    output_file_kilosort = [output_dir,'kilosort_',num2str(tag),'_probe',num2str(probe),'shank',num2str(shank),'_data.bin'];
    fidout = fopen(output_file_kilosort, 'w'); 
    
    TrialNum_prev = nan;
    TrialNum = nan;
    Trial_time_all=[];
    for i_trial = 1:size(input_filelist,1)
        disp('--------------------------------------------');
        disp(['processing file ',input_filelist(i_trial).name]);
        read_Intan_RHD2000_file_auto(input_filelist(i_trial).name);
        
        ch_all_raw = amplifier_data';  % amplifier_data is units of microvolts
        dig_all_raw = board_dig_in_data';
    
        VoltageTraceInV_allCh = ch_all_raw;
        clear ch_all_raw
        
        Trigger_allCh = dig_all_raw(:,trigger_ch);
        if bitcode_ch ~= 0
            Bitcode_allCh = dig_all_raw(:,bitcode_ch);
        else
            Bitcode_allCh = [];
        end
        clear ch_all_raw dig_all_raw
        
        TimeStamps = t_amplifier;
        idx=find(Bitcode_allCh>0.2);
        if TimeStamps(idx(1))>1000
            t_sample = 1.57;
        else
            t_sample = 0.57;
        end
    
        ch_MUA=[]; 
        if probe==1
            ch_MUA = VoltageTraceInV_allCh(:,1:64);
        else
            ch_MUA = VoltageTraceInV_allCh(:,65:128);
        end
        clear VoltageTraceInV_allCh
        
        for i_ch = 1:ch_spike
            ch_tmp         = timeseries(ch_MUA(:,i_ch),TimeStamps);
            ch_tmp_MUA     = idealfilter(ch_tmp,[300 6000],'pass');
            ch_MUA(:,i_ch) = ch_tmp_MUA.data;
        end
        
        commonNoise = trimmean(ch_MUA,40,2);
        i_noise=[find(commonNoise>150);find(commonNoise<-200)]; 
        for i_ch=1:ch_spike
            t_post_stim = [0 t_sample+1.3]; 
            i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
            idx1=i_post_stim(end);
            X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
            b1 = regress(ch_MUA(i_post_stim,i_ch),X);
            t_post_stim = [t_sample+1.3 t_sample+2.1];
            i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
            idx2=i_post_stim(end);
            X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
            b2 = regress(ch_MUA(i_post_stim,i_ch),X);
            t_post_stim = [t_sample+2.1 t_sample+4.5]; 
            i_post_stim = find(TimeStamps>t_post_stim(1) & TimeStamps<t_post_stim(2));
            X = [ones(size(commonNoise(i_post_stim),1),1) commonNoise(i_post_stim)];
            b3 = regress(ch_MUA(i_post_stim,i_ch),X);
            ch_MUA(:,i_ch) = ch_MUA(:,i_ch) - [commonNoise(1:idx1)*b1(2);commonNoise(idx1+1:idx2)*b2(2);commonNoise(idx2+1:end)*b3(2)];
        end

        if ~isempty(i_noise)
            for j=1:length(i_noise)
                if i_noise(j)>1&&i_noise(j)<length(commonNoise)-1
                    ch_MUA(i_noise(j)-1:i_noise(j)+1,:)=zeros(3,size(ch_MUA,2));
                end
            end
        end
        clear commonNoise
        
        % sort channels
        ch_MUA = func_sortChannel_DBC64ch(ch_MUA);  
        if shank==1
           ch_MUA = ch_MUA(:,1:32);
        else
           ch_MUA = ch_MUA(:,33:64);
        end

        % read bit code
        if max(Bitcode_allCh(2:end-1))>0.2
            if bitcode_ch ~= 0&isempty(strfind(input_filelist(i_trial).name,'passive'))
                TrialNum_prev = TrialNum;
                try
                    TrialNum = func_read_bitcode(Bitcode_allCh,TimeStamps);
                    disp(['done!  Matched to solo trial #',num2str(TrialNum)]);
                catch
                    TrialNum = TrialNum_prev+1;
                    warning(['Bitcode Failed!  Assigned to solo trial #',num2str(TrialNum)]);
                    keyboard
                end
            elseif bitcode_ch ~= 0&~isempty(strfind(input_filelist(i_trial).name,'passive'))
                TrialNum=-func_read_bitcode(Bitcode_allCh,TimeStamps);
            else
                TrialNum = -i_trial;
            end
        
        else
            TrialNum=0;
        end
    
        Trial_time_tmp=[TimeStamps;ones(1,length(TimeStamps))*TrialNum];
        Trial_time_all=[Trial_time_all Trial_time_tmp];
        
        if ismember(TrialNum,laser_trial)
             output_file_name = [output_dir,'raw_trace_',num2str(tag),'_probe',num2str(probe),'shank',num2str(shank),'_laser_trial_',num2str(TrialNum),'.mat'];
             disp(['saving laser trial: ',output_file_name]);
             save(output_file_name,'ch_MUA','TimeStamps','Trigger_allCh','Bitcode_allCh');
        end
        
        disp(['write: trial ',num2str(TrialNum)]);
        ch_MUA = int16(ch_MUA');
        fwrite(fidout, ch_MUA, 'int16');
        clear ch_MUA ch_tmp_MUA VoltageTraceInV_allCh TimeStamps Trial_time_tmp Trigger_allCh allOther_allCh Bitcode_allCh
    end 
    fclose(fidout);
    
    output_file=[output_dir,'Trial_time_all_',num2str(tag),'_probe',num2str(probe),'shank',num2str(shank),'_trial_',num2str(Trial_time_all(2,1)),'_',num2str(Trial_time_all(2,end)),'.mat'];
    save(output_file,'Trial_time_all','-v7.3');
    events=find(Trial_time_all(1,:)==0)'/20000;
    csvwrite([output_dir,'events.csv'],events);
return
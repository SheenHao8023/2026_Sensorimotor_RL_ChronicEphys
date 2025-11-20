%% 原始 Intan 数据 → Kilosort spike sorting → 提取刺激事件
clear
datapath={
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_06\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_07\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_08\Intan\';...
    'C:\Users\XinHao\Desktop\BAYLORGC130\2021_02_09\Intan\';...
    };

Laser_trial={
    [3 4 5];...
    [3 4 5];...
    [3 4 5];...
    [3 4 5];...
    };%sti right/left/bilat

for k=1:length(datapath)
    laser_trial=Laser_trial{k};      %sti right/left/bilat
    file_dir = datapath{k};
    file_list_tmp = dir([file_dir,'/*.rhd']);   % 找到目标文件
    file_list=[];
    for j=1:size(file_list_tmp,1)  % for passive DBC64ch recording
        if ~isempty(strfind(file_list_tmp(j).name,'GC130'))  %寻找包含命名的文件
            file_list(end+1).name=file_list_tmp(j).name;
        end
    end

    input_filelist = [];
    for i_file = 3:size(file_list,2)-1  %排除前两个和最后一个文件，photo stim trial to check noise
        input_filelist(end+1).name = [file_dir,file_list(i_file).name];
    end

    for probe=1:2 %left ALM 1 intan port A 1-64 ch;right ALM 2 intan port B 65-128 ch
        for shank=1:2 %1 2
            if ~exist([file_dir,'Kilosort_data',num2str(probe),num2str(shank)] ,'dir')
                mkdir([file_dir,'Kilosort_data',num2str(probe),num2str(shank)])
            end
            voltage_dir = [file_dir,'Kilosort_data',num2str(probe),num2str(shank),'\']; %作为kilosort输入
            func_process_voltage_traces_kilosort(input_filelist, voltage_dir, 64, 1, 2, probe, shank, laser_trial);

            folder=[file_dir,'Kilosort_data',num2str(probe),num2str(shank),'\']; 
            matfile=dir([folder,'Trial_time*.mat']);
            load([folder,matfile(1).name]);
            events=find(Trial_time_all(1,:)==0)';
            writematrix(events, fullfile(folder, 'events.csv'));
            clear Trial_time_all;

        end
    end
    % master_kilosort;   % use Kilosort 2.0
end
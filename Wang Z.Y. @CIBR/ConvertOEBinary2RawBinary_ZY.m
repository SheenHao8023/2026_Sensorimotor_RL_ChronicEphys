clear
addpath(genpath('C:\Kilosort2-user\my_kilosort'));
datapath={
    'C:\Users\XinHao\Desktop\CIBRZY25\2025-10-12_15-40-45\';...
    };

Laser_trial={
    [1 2 3];...
    [1 2 3];...
    [1 2 3];...
    [1 2 3];...
    };%sti right/left/bilat

for k=1:length(datapath)
    laser_trial=Laser_trial{k};%sti right/left/bilat
    file_dir = datapath{k};
    i_str=strfind(file_dir,'\');

    % for multiple neuropixels probe
    searchPath = [file_dir ,'**\settings.xml']; % Search in folder and subfolders for  *.xml
    Files      = dir(searchPath); % Find all .xml files
    XMLfile=[Files(1).folder,'\',Files(1).name];%
    searchPath = [file_dir ,'**\experiment1\recording1\structure.oebin']; % Search in folder and subfolders for  *.oebin
    Files      = dir(searchPath); % Find all .oebin files
    nprobe=length(Files);
    npType=[];
    probek=0;% default: run probe one by one if you need, probek=0 run all probes together; probek=1 2 3... run the "probek" th probe.

    mmapN=1;
    oepxi=1;

    % no signal
    if mmapN==1
        probe=nprobe;%for multiple NP
        shank = 1; %CB H2 probe 32*2 shank1 1-32, shank2 33-64
        for probe_tmp=1:probe
            for shank_tmp=1:shank% NP1 only 1 shank, NP2 4 shank treat as 1 shank
                bin_path=[file_dir,'Kilosort_data',num2str(probe_tmp),num2str(shank_tmp)];
                if ~exist(bin_path)%
                    mkdir(bin_path)%
                end
                [xcoords, ycoords, nptype]=createChannelMapFile_384_NPx(XMLfile,probe_tmp,bin_path);
                npType(probe_tmp)=nptype;%for multiple NP1 or 2 recording
            end
        end
        func_process_voltage_traces_kilosort_OE_NPx(file_dir, 384, probe, shank, laser_trial,oepxi, probek);%for multiple NP1 or 2 recording
    else
        probe=1; %region 1 intan headstage 1-64 ch;region 2 intan headstage 65-128 ch etc
        shank=2;%CB H2 probe 32*2 shank1 1-32, shank2 33-64
        % generate output file directorys
        for probe_tmp=1:probe
            for shank_tmp=1:shank
                if ~exist([file_dir,'Kilosort_dataOE',num2str(probe_tmp),num2str(shank_tmp)])%
                    mkdir([file_dir,'Kilosort_dataOE',num2str(probe_tmp),num2str(shank_tmp)])%
                end
            end
        end
        func_process_voltage_traces_kilosort_OE_CB64(file_dir, 32, probe, shank, laser_trial, mmapN);
    end
    % master_kilosort_NPx_ZY; 
end

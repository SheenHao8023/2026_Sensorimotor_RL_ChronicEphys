clear
datapath={
    'D:\Data\ZYW\2025-10-12_15-40-45\';...
    };
Laser_trial={
    [1 2 3];...
    [1 2 3];...
    [1 2 3];...
    [1 2 3];...
    };%sti right/left/bilat

for k=1:length(datapath)

    % for multiple neuropixels probe
    file_dir = datapath{k};
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
    if mmapN==1
        probe=nprobe; %for multiple NP
    else
        probe=1; %region 1 intan headstage 1-64 ch;region 2 intan headstage 65-128 ch etc
    end

    if probek==0
        runprobe=1:probe;
        shank = 1; %CB H2 probe 32*2 shank1 1-32, shank2 33-64
    else
        runprobe=probek;
        shank = 2; %CB H2 probe 32*2 shank1 1-32, shank2 33-64
    end

    for probe_tmp=runprobe
        for shank_tmp=1:shank
            if mmapN==1
                rootZ = [datapath{k},'Kilosort_data',num2str(probe_tmp),num2str(shank_tmp),'\']; %,num2str(position) the raw data binary file is in this folder
            else
                rootZ = [datapath{k},'Kilosort_dataOE',num2str(probe_tmp),num2str(shank_tmp),'\']; %,num2str(position) the raw data binary file is in this folder
            end
            addpath(genpath('D:\Kilosort2-master')) % path to kilosort folder
            addpath('D:\Kilosort2-user\my_kilosort\npy-matlab-master\npy-matlab') % for converting to Phy
            rootH = rootZ; % path to temporary binary file (same size as data, should be on fast SSD)
            pathToYourConfigFile = 'D:\Kilosort2-user\my_kilosort\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
            [~, ~, nptype]=createChannelMapFile_384_NPx(XMLfile, probe_tmp, rootZ);
            npType(probe_tmp) = nptype;
            if npType(probe_tmp)==1
                chanMapFile = 'neuropixPhase3A_kilosortChanMap_NP1.mat';
            else
                chanMapFile = 'neuropixels_kilosortChanMap_NP2_4shank.mat';
            end

            ops.trange = [0 Inf]; % time range to sort
            ops.NchanTOT    = 384; % total number of channels in your recording

            if npType(probe_tmp)==1
                run(fullfile(pathToYourConfigFile, 'configFile384_NP1.m'))
            else
                run(fullfile(pathToYourConfigFile, 'configFile384_NP2.m'))
            end
            ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
            ops.chanMap = fullfile(rootZ, chanMapFile);

            %%% this block runs all the steps of the algorithm
            fprintf('Looking for data inside %s \n', rootZ)

            % is there a channel map file in this folder?
            fs = dir(fullfile(rootZ, 'chan*.mat'));
            if ~isempty(fs)
                ops.chanMap = fullfile(rootZ, fs(1).name);
            end

            % find the binary file
            fs  = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
            ops.fbinary = fullfile(rootZ, fs(1).name);

            % preprocess data to create temp_wh.dat
            rez = preprocessDataSub(ops);

            % time-reordering as a function of drift
            rez = clusterSingleBatches(rez);

            % saving here is a good idea, because the rest can be resumed after loading rez
            save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

            % main tracking and template matching algorithm
            rez = learnAndSolve8b(rez);

            % final merges
            rez = find_merges(rez, 1);

            % final splits by SVD
            rez = splitAllClusters(rez, 1);

            % final splits by amplitudes
            rez = splitAllClusters(rez, 0);

            % decide on cutoff
            rez = set_cutoff(rez);

            fprintf('found %d good units \n', sum(rez.good>0))

            % write to Phy
            fprintf('Saving results to Phy  \n')
            rezToPhy(rez, rootZ);

            clearvars -except datapath Laser_trial k position probe shank mmapN npType probek runprobe XMLfile searchPath Files  

        end
    end
end
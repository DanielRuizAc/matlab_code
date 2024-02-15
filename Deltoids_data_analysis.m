% dinfo = dir('./OriginalData/DeltBicep/NoExo/*.txt');
dinfo_without = dir('./OriginalData/FnDeltBiceps/withoutExo/*.txt');
load(".\DataSets\OnlyXsens\kmeans_data.mat")

without_exo_data = dictionary;
Desired_segments = {'LeftForeArm', 'RightForeArm'};
Desired_joints = {'RightShoulder_RightUpperArm','LeftShoulder_LeftUpperArm'};

transform = true;
% transform = false;

fr=100;
nq_freq = fr/2;
cutoff_low = 3.0;

[b,a] = butter(4, cutoff_low/nq_freq, "low");

for ind_file = 1 : length(dinfo_without)
    fileName = strcat('./OriginalData/FnDeltBiceps/withoutExo/', dinfo_without(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, true);
    discarted_frames

    exp_data.original_frames = exp_frames;
    modified_frames = CreateDBXSens(exp_frames,Desired_segments, Desired_joints, transform);

    desired_fields.position=[1 3];
    Desired_vars_joints = [1 2 3];

    unified_db = UnifyCharsXSens(modified_frames, desired_fields, Desired_vars_joints);
    dat_filt = unified_db;
    dims = size(dat_filt);
    for i_fil=1:dims(2)
        dat_filt(:,i_fil) = filtfilt(b,a,dat_filt(:,i_fil));
    end

    % exp_data.filtered_db = dat_filt;
    norm_db = (dat_filt - C).*(1./S);

    for i=1:size(norm_db,1)
        dataset_pca(i,:) = coeff_pca'*(norm_db(i,:)');
    end
    exp_data.pca_db = dataset_pca;
    without_exo_data(dinfo_without(ind_file).name) = exp_data;
end


dinfo_with = dir('./OriginalData/FnDeltBiceps/withExo/*.txt');
with_exo_data = dictionary;
for ind_file = 1 : length(dinfo_with)
    fileName = strcat('./OriginalData/FnDeltBiceps/withExo/', dinfo_with(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, true);
    discarted_frames

    exp_data.original_frames = exp_frames;
    modified_frames = CreateDBXSens(exp_frames,Desired_segments, Desired_joints, transform);

    desired_fields.position=[1 3];
    Desired_vars_joints = [1 2 3];

    unified_db = UnifyCharsXSens(modified_frames, desired_fields, Desired_vars_joints);
    dat_filt = unified_db;
    dims = size(dat_filt);
    for i_fil=1:dims(2)
        dat_filt(:,i_fil) = filtfilt(b,a,dat_filt(:,i_fil));
    end

    % exp_data.filtered_db = dat_filt;
    norm_db = (dat_filt - C).*(1./S);

    for i=1:size(norm_db,1)
        dataset_pca(i,:) = coeff_pca'*(norm_db(i,:)');
    end
    exp_data.pca_db = dataset_pca;
    with_exo_data(dinfo_with(ind_file).name) = exp_data;
end

% save("DataSets\withoutExoInertial", "without_exo_data")
% save("DataSets\withExoInertial", "with_exo_data")

%% Clustering
dinfo_with = dir('./OriginalData/FnDeltBiceps/withExo/*.txt');
dinfo_without = dir('./OriginalData/FnDeltBiceps/withoutExo/*.txt');

clustering_without_exo = dictionary;
for ind_file = 1 : length(dinfo_without)
    exp_data = without_exo_data(dinfo_without(ind_file).name);
    [~,idx_test] = pdist2(C_eucl,exp_data.pca_db(:,1),'euclidean','Smallest',1);
    clust_data = struct;
    clust_data.ids = idx_test;
    clust_data.pca = exp_data.pca_db(:,1); 
    clustering_without_exo(dinfo_without(ind_file).name) = clust_data;
end

clustering_with_exo = dictionary;
for ind_file = 1 : length(dinfo_with)
    exp_data = with_exo_data(dinfo_with(ind_file).name);
    [~,idx_test] = pdist2(C_eucl,exp_data.pca_db(:,1),'euclidean','Smallest',1);
    clust_data = struct;
    clust_data.ids = idx_test;
    clust_data.pca = exp_data.pca_db(:,1);
    clustering_with_exo(dinfo_with(ind_file).name) = clust_data;
end

save("DataSets\clustering_data", "clustering_with_exo", "clustering_without_exo")

%%

dinfo_with = dir('./OriginalData/FnDeltBiceps/withExo/*.txt');
dinfo_without = dir('./OriginalData/FnDeltBiceps/withoutExo/*.txt');
fr = 1000;
nq_freq = fr/2;
cutoff_low = 10;
cutoff_high = 200;

[b,a] = butter(4,[cutoff_low,cutoff_high]/nq_freq,"bandpass");

max_emg_vals = struct;
emg_without_exo_data = dictionary;
for ind_file = 1 : length(dinfo_without)
    emg_exp_data = struct;
    fileName = strcat('./OriginalData/FnDeltBiceps/withoutExo/', dinfo_without(ind_file).name);
    % if (without_exo_data.isKey(fileName))
    %     continue
    % end
    [channels_data] = ObtainEMGData(fileName, true);
    fieldnames_emg = fieldnames(channels_data);
    emg_exp_data.ch_data = channels_data;
    for i=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{i};
        if (isfield(max_emg_vals, ch_name)==false)
            max_emg_vals.(ch_name) = 0;
        end
        channels_data.(ch_name)(2,:) = abs(filtfilt(b,a,channels_data.(ch_name)(2,:)));
        max_emg_vals.(ch_name) = max(max_emg_vals.(ch_name), max(channels_data.(ch_name)(2,:)));
    end
    % channels_data(2,:) = abs(filtfilt(b,a,channels_data(2,:)));
    emg_exp_data.rect_filt = channels_data;
    emg_without_exo_data(dinfo_without(ind_file).name) = emg_exp_data;
end


emg_with_exo_data = dictionary;
for ind_file = 1 : length(dinfo_with)
    emg_exp_data = struct;
    fileName = strcat('./OriginalData/FnDeltBiceps/withExo/', dinfo_with(ind_file).name);
    % if (without_exo_data.isKey(fileName))
    %     continue
    % end
    [channels_data] = ObtainEMGData(fileName, true);
    fieldnames_emg = fieldnames(channels_data);
    emg_exp_data.ch_data = channels_data;
    for i=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{i};
        if (isfield(max_emg_vals, ch_name))
            max_emg_vals.(ch_name) = 0;
        end
        channels_data.(ch_name)(2,:) = abs(filtfilt(b,a,channels_data.(ch_name)(2,:)));
        max_emg_vals.(ch_name) = max(max_emg_vals.(ch_name), max(channels_data.(ch_name)(2,:)));
    end
    % channels_data(2,:) = abs(filtfilt(b,a,channels_data(2,:)));
    emg_exp_data.rect_filt = channels_data;
    emg_with_exo_data(dinfo_with(ind_file).name) = emg_exp_data;
end

save("DataSets\withoutExoEMG", "emg_without_exo_data")
save("DataSets\withExoEMG", "emg_with_exo_data")
save("DataSets\maxEMGVals", "max_emg_vals")
%%
% fieldnames_exps = fieldnames(emg_exp_data);
% for i=1:length(fieldnames_exps)
%     fld_name = fieldnames_exps{i};
%     emg_without_exo_data.(fld_name).rect_filt = 
% end
% 
% 
% %%
% exp_data = without_exo_data("2023-11-22-133039.txt");
% emg_exp_data = emg_without_exo_data("2023-11-22-133039.txt");
% 
% [~,idx_test] = pdist2(C_eucl,exp_data.pca_db(:,1),'euclidean','Smallest',1);
% valid_indexes = find((idx_test==2));
% recorted_data = struct;
% 
% for i=1:length(fieldnames_emg)
%     ch_name = fieldnames_emg{i};
%     channel_data = emg_exp_data.ch_data.(ch_name);
%     emg_indexes = round(channel_data(1,:)/10);
%     recorted_data.(ch_name).recorted = channel_data(2, ismember(emg_indexes, valid_indexes));
%     WinLen = 2000;                                            % Window Length For RMS Calculation
%     index=1;
%     step=200;
%     window_indexes = 1:step:length(recorted_data.(ch_name).recorted)-WinLen;
% 
%     rmsv = zeros(1, length(window_indexes));
%     arvv = zeros(1, length(window_indexes));
%     md_freq = zeros(1, length(window_indexes));
%     mean_freq = zeros(1, length(window_indexes));
%     index=1;
%     for i=window_indexes
%         signal_cutted = recorted_data.(ch_name).recorted(i:i+WinLen);
%         rmsv(index) = rms(signal_cutted,2);
%         arvv(index) = sum(abs(signal_cutted),2)/WinLen;
%         md_freq(index) = medfreq(signal_cutted,fr);
%         mean_freq(index) = meanfreq(signal_cutted,fr);
%         index=index+1;
%     end
% 
%     recorted_data.(ch_name).rmsv = rmsv;
%     recorted_data.(ch_name).arvv = arvv;
%     recorted_data.(ch_name).md_freq = md_freq;
%     recorted_data.(ch_name).mean_freq = mean_freq;
% end


% emg_idexes = round(emg_exp_data.ch_data.ch_0(1,:)/10);
% valid_emg_indexes = find(ismember(emg_idexes, valid_indexes));
% 
% EMG_valid_indexes = emg_exp_data.ch_data.ch_0(1,(emg_exp_data.ch_data.ch_0(1,:)));




% WinLen = 2000;                                            % Window Length For RMS Calculation
% index=1;
% step=2000;
% 
% without_exo_keys = without_exo_data.keys;
% for key_id=1:length(without_exo_keys)
%     key = without_exo_keys(key_id);
%     filt_ch_data = without_exo_data(key).ch_data; % without_exo_data(key).rect_filt;
% 
%     ch_names = fieldnames(filt_ch_data);
% 
%     fat_data = struct;
%     for ch_ind=1:length(ch_names)
%         ch_nm = ch_names{ch_ind};
%         ch_data = filt_ch_data.(ch_nm);
%         rmsv = [];
%         md_freq = [];
%         mean_freq = [];
%         index = 1;
%         for i=1:step:size(ch_data,2)-WinLen
%             signal_cutted= ch_data(2,i:i+WinLen);
%             rmsv(index) = rms(signal_cutted,2);
%             md_freq(index) = medfreq(signal_cutted,fr);
%             mean_freq(index) = meanfreq(signal_cutted,fr);
%             index=index+1;
%         end
%         fat_data.(ch_nm).rms_dat = rmsv;
%         fat_data.(ch_nm).md_freq = md_freq;
%         fat_data.(ch_nm).mn_freq = mean_freq;
%     end
%     without_exo_data(key).fat_data = fat_data;
% end

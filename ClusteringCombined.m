dinfo = dir('./OriginalData/XSensEmg/WithoutExo/*.txt');
load(".\DataSets\OnlyXsens\kmeans_data.mat")
without_exo_data = dictionary;

Desired_segments = {'LeftForeArm', 'RightForeArm'};
Desired_joints = {'RightShoulder_RightUpperArm','LeftShoulder_LeftUpperArm'};

transform = true;
% transform = false;

fr=100;
nq_freq = fr/2;
cutoff_low = 3.0;

[b,a] = butter(4, cutoff_low/nq_freq,"low");


for ind_file = 1 : length(dinfo)
    fileName = strcat('./OriginalData/XSensEmg/WithoutExo/', dinfo(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, true);
    discarted_frames

    exp_data.original_frames = exp_frames;
    exp_data.modified_frames = CreateDBXSens(exp_frames,Desired_segments, Desired_joints, transform);

    desired_fields.position=[1 3];
    Desired_vars_joints = [1 2 3];

    exp_data.unified_db = UnifyCharsXSens(exp_data.modified_frames, desired_fields, Desired_vars_joints);
    dat_filt = exp_data.unified_db;
    dims = size(dat_filt);
    for i_fil=1:dims(2)
        dat_filt(:,i_fil) = filtfilt(b,a,dat_filt(:,i_fil));
    end

    exp_data.filtered_db = dat_filt;
    exp_data.norm_db = (dat_filt - C).*(1./S);

    for i=1:size(exp_data.norm_db,1)
        dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
    end
    exp_data.pca_db = dataset_pca;

    without_exo_data(fileName) = exp_data;
end
%%

dinfo = dir('./OriginalData/XSensEmg/WithExo/*.txt');
with_exo_data = dictionary;

for ind_file = 1 : length(dinfo)
    fileName = strcat('./OriginalData/XSensEmg/WithExo/', dinfo(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, true);
    discarted_frames
    exp_data.original_frames = exp_frames;
    exp_data.modified_frames = CreateDBXSens(exp_frames,Desired_segments, Desired_joints, transform);

    desired_fields.position=[1 3];
    Desired_vars_joints = [1 2 3];

    exp_data.unified_db = UnifyCharsXSens(exp_data.modified_frames, desired_fields, Desired_vars_joints);
    
    dat_filt = exp_data.unified_db;
    dims = size(dat_filt);
    for i_fil=1:dims(2)
        dat_filt(:,i_fil) = filtfilt(b,a,dat_filt(:,i_fil));
    end

    exp_data.filtered_db = dat_filt;
    exp_data.norm_db = (exp_data.filtered_db - C).*(1./S);
    
    for i=1:size(exp_data.norm_db,1)
        dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
    end    
    exp_data.pca_db = dataset_pca;
    
    with_exo_data(fileName) = exp_data;
end


%% Clustering

without_exo_data = without_exo_data.remove("./OriginalData/XSensEmg/WithoutExo/2023-10-27-103755.txt");

with_exo_data = with_exo_data.remove("./OriginalData/XSensEmg/WithExo/2023-10-27-103755.txt");

without_exo_keys = without_exo_data.keys;


without_exo_clustered = dictionary;

for key_id=1:length(without_exo_keys)
    key = without_exo_keys(key_id);
    exp_data = without_exo_data(key);
    [~,idx_test] = pdist2(C_eucl,exp_data.pca_db(:,1),'euclidean','Smallest',1);
    % data_clust = zeros(length(idx_test),2);
    data_clust = struct;
    data_clust.frame_number = exp_data.original_frames.frame_number;
    data_clust.labels = idx_test;
    % figure
    % plot(idx_test)

    without_exo_clustered(key) = data_clust;
end

with_exo_keys = with_exo_data.keys;
with_exo_clustered = dictionary;

for key_id=1:length(with_exo_keys)
    key = with_exo_keys(key_id);
    exp_data = with_exo_data(key);
    [~,idx_test] = pdist2(C_eucl,exp_data.pca_db(:,1),'euclidean','Smallest',1);
    % figure
    % plot(idx_test)
    data_clust = struct;
    data_clust.frame_number = exp_data.original_frames.frame_number;
    data_clust.labels = idx_test;

    with_exo_clustered(key) = data_clust;
end


%% Load EMG data

load("DataSets\OnlyEMG\EMG_DATA.mat")


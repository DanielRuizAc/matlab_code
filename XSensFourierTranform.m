load("DataSets\OnlyXsens\transTrainData.mat");

n_exps_no_exo = length(without_exo_data.keys);
n_exps_exo = length(with_exo_data.keys);
n_max_exps = max(n_exps_exo,n_exps_no_exo);
without_exo_keys = without_exo_data.keys;

fr = 100;

exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
transf_data = exp_data.modified_frames;

r_forearm_data = transf_data.segments.RightForeArm;

dims_fields = fieldnames(r_forearm_data);
for i=1:length(dims_fields)
    fld_name = dims_fields(i);
    data=r_forearm_data.(fld_name{1})(:,1);
    f = fr*(0:length(data)/2)/length(data);
    four = fft(data);
    four = abs(four/length(four));
    four = four(1:(length(four)/2+1));
    four(2:end-1) = 2*four(2:end-1);
    if (fld_name=="ang_vel")
        fld_name="angular velocity";
    elseif (fld_name=="ang_accel")
        fld_name="angular acceleration";
    end
    figure
    plot(f, four)

    title(strcat("Fourier transform ",fld_name," ForeArm"))
    xlabel("Frequency Hz")
    ylabel("|fft|")
end

shoulder_data = transf_data.joints.RightShoulder_RightUpperArm;
data=shoulder_data(:,1);
f = fr*(0:length(data)/2)/length(data);
four = fft(data);
four = abs(four/length(four));
four = four(1:(length(four)/2+1));
four(2:end-1) = 2*four(2:end-1);

figure
plot(f, four)
title(strcat("Fourier transform right shoulder joint"))
xlabel("Frequency Hz")
ylabel("|fft|")


%% 
load("DataSets\OnlyXsens\transTrainData.mat");

exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
transf_data = exp_data.modified_frames;

r_forearm_data = transf_data.segments.RightForeArm;

figure(1)
subplot(2,2,2)
plot((1:length(r_forearm_data.position(:,1)))*0.01, r_forearm_data.position(:,1))
title(strcat("x component of position with respect to pelvis segment"))
xlabel("time (s)")
ylabel("position (m)")

subplot(2,2,4)
plot((1:length(r_forearm_data.position(:,2)))*0.01, r_forearm_data.position(:,2))
title(strcat("y component of position with respect to pelvis segment"))
xlabel("time (s)")
ylabel("position (m)")

load("DataSets\OnlyXsens\noTransTrainData.mat");

exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
transf_data = exp_data.modified_frames;

r_forearm_data = transf_data.segments.RightForeArm;

figure(1)
subplot(2,2,1)
plot((1:length(r_forearm_data.position(:,1)))*0.01, r_forearm_data.position(:,1))
title(strcat("x component of position with respect to base frame"))
xlabel("time (s)")
ylabel("position (m)")

subplot(2,2,3)
plot((1:length(r_forearm_data.position(:,2)))*0.01, r_forearm_data.position(:,2))
title(strcat("y component of position with respect to base frame"))
xlabel("time (s)")
ylabel("position (m)")


%% 
load("DataSets\OnlyXsens\transTrainData.mat");
load("DataSets\OnlyXsens\kmeans_data.mat")
exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
dat_filt = exp_data.filtered_db;

exp_data.norm_db = (dat_filt - C).*(1./S);

figure
for i=1:10
hold on
plot((1:length(exp_data.norm_db(:,i)))*0.01, exp_data.norm_db(:,i), "LineWidth", 1.2)
end
title("Normalized characteristics")
xlabel("time (s)")
ylabel("characteristic normalized amplitude")

legend(["Right forearm x", "Right forearm z", "Left forearm x", "Left forearm z", ...
    "Right Shoulder rotation x", "Right Shoulder rotation y", "Right Shoulder rotation z", ...
    "Left Shoulder rotation x", "Left Shoulder rotation y", "Left Shoulder rotation z",])


%%
clear dataset_pca

load("DataSets\OnlyXsens\transTrainData.mat");
load("DataSets\OnlyXsens\kmeans_data.mat")
exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
dat_filt = exp_data.filtered_db;

exp_data.norm_db = (dat_filt - C).*(1./S);
for i=1:size(exp_data.norm_db,1)
    dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
end
exp_data.pca_db = dataset_pca;

figure
subplot(3,1,1)
plot((1:length(dataset_pca(:,1)))*0.01, dataset_pca(:,1), "LineWidth", 1.2)
title("Transformed Characteristic 1")
xlabel("time (s)")
ylabel("Value")
subplot(3,1,2)
plot((1:length(dataset_pca(:,2)))*0.01, dataset_pca(:,2), "LineWidth", 1.2)
title("Transformed Characteristic 2")
xlabel("time (s)")
ylabel("Value")
subplot(3,1,3)
plot((1:length(dataset_pca(:,3)))*0.01, dataset_pca(:,3), "LineWidth", 1.2)
title("Transformed Characteristic 3")
xlabel("time (s)")
ylabel("Value")


exp_data = with_exo_data("./OriginalData/OnlyXsens/WithExo/20231017 150656.271.txt");
dat_filt = exp_data.filtered_db;

exp_data.norm_db = (dat_filt - C).*(1./S);
for i=1:size(exp_data.norm_db,1)
    dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
end
exp_data.pca_db = dataset_pca;

subplot(3,1,1)
hold on
plot((1:length(dataset_pca(:,1)))*0.01, dataset_pca(:,1), "LineWidth", 1.2)
legend(["without exoskeleton", "with exoskeleton"])
subplot(3,1,2)
hold on
plot((1:length(dataset_pca(:,2)))*0.01, dataset_pca(:,2), "LineWidth", 1.2)
legend(["without exoskeleton", "with exoskeleton"])
subplot(3,1,3)
hold on
plot((1:length(dataset_pca(:,3)))*0.01, dataset_pca(:,3), "LineWidth", 1.2)
legend(["without exoskeleton", "with exoskeleton"])


%% 


load("DataSets\OnlyXsens\transTrainData.mat");
load("DataSets\OnlyXsens\kmeans_data.mat")
exp_data = without_exo_data("./OriginalData/OnlyXsens/WithoutExo/20231017 144722.296.txt");
dat_filt = exp_data.filtered_db;

exp_data.norm_db = (dat_filt - C).*(1./S);
for i=1:size(exp_data.norm_db,1)
    dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
end
exp_data.pca_db = dataset_pca;

[~,idx_test] = pdist2(C_eucl, dataset_pca(:,1),'euclidean','Smallest',1);
rest_indexes = idx_test==1;
task_indexes = idx_test==2;

[~ ,from_rest_to_work] = max(diff(idx_test));
[~ ,from_work_to_rest] = min(diff(idx_test));


t=(1:length(dataset_pca(:,1)))*0.01;
fig = figure;
subplot(2,1,1)
plot(t, dataset_pca(:,1), "LineWidth", 1.2)
hold on
plot(t, ones(1,length(dataset_pca(:,1))).*C_eucl(1), "LineWidth", 1.0)
hold on
plot(t, ones(1,length(dataset_pca(:,1))).*C_eucl(2), "LineWidth", 1.0)
ylim([min(dataset_pca(:,1)) - 1.3, max(dataset_pca(:,1)) + 0.7])
hold on
xline(from_rest_to_work*0.01,'r-',{'From rest to', 'task transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
hold on
xline(from_work_to_rest*0.01,'r-',{'From task to', 'rest transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);

hold on
xlabel("time (s)")
ylabel("Characteristic Value")
legend(["PCA Characteristic", "'rest' centroid value", "'task' centroid value"])
title("Experiment without exoskeleton")
subplot(2,1,2)
plot(t(rest_indexes), idx_test(rest_indexes), 'o')
hold all
plot(t(task_indexes), idx_test(task_indexes), 'o')
ylim([0.7, 2.3])
xlabel("time (s)")
ylabel("State")
hold on
xline(from_rest_to_work*0.01,'r-',{'From rest to', 'task transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
hold on
xline(from_work_to_rest*0.01,'r-',{'From task to', 'rest transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
legend(["'Rest' state", "'Task' state"])


exp_data = with_exo_data("./OriginalData/OnlyXsens/WithExo/20231017 150656.271.txt");
dat_filt = exp_data.filtered_db;

exp_data.norm_db = (dat_filt - C).*(1./S);
for i=1:size(exp_data.norm_db,1)
    dataset_pca(i,:) = coeff_pca'*(exp_data.norm_db(i,:)');
end
exp_data.pca_db = dataset_pca;

[~,idx_test] = pdist2(C_eucl, dataset_pca(:,1),'euclidean','Smallest',1);
rest_indexes = idx_test==1;
task_indexes = idx_test==2;

[~ ,from_rest_to_work] = max(diff(idx_test));
[~ ,from_work_to_rest] = min(diff(idx_test));


t=(1:length(dataset_pca(:,1)))*0.01;
figure
subplot(2,1,1)
plot(t, dataset_pca(:,1), "LineWidth", 1.2)
hold on
plot(t, ones(1,length(dataset_pca(:,1))).*C_eucl(1), "LineWidth", 1.0)
hold on
plot(t, ones(1,length(dataset_pca(:,1))).*C_eucl(2), "LineWidth", 1.0)
ylim([min(dataset_pca(:,1)) - 1.3, max(dataset_pca(:,1)) + 0.7])
hold on
xline(from_rest_to_work*0.01,'r-',{'From rest to', 'task transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
hold on
xline(from_work_to_rest*0.01,'r-',{'From task to', 'rest transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);

hold on
xlabel("time (s)")
ylabel("Characteristic Value")
legend(["PCA Characteristic", "'rest' centroid value", "'task' centroid value"])
title("Experiment with exoskeleton")
subplot(2,1,2)
plot(t(rest_indexes), idx_test(rest_indexes), 'o')
hold on
plot(t(task_indexes), idx_test(task_indexes), 'o')
ylim([0.7, 2.3])
xlabel("time (s)")
ylabel("State")
hold on
xline(from_rest_to_work*0.01,'r-',{'From rest to', 'task transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
hold on
xline(from_work_to_rest*0.01,'r-',{'From task to', 'rest transition'}, ...
    "LabelOrientation", "horizontal", ...
    "LabelVerticalAlignment", "bottom", ...
    "LineWidth", 0.9);
legend(["'Rest' state", "'Task' state"])
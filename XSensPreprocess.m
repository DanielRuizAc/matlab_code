clear

dinfo = dir('./OriginalData/OnlyXsens/WithoutExo/*.txt');
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
    fileName = strcat('./OriginalData/OnlyXsens/WithoutExo/', dinfo(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, false);
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
    without_exo_data(fileName) = exp_data;
end


dinfo = dir('./OriginalData/OnlyXsens/WithExo/*.txt');
with_exo_data = dictionary;

for ind_file = 1 : length(dinfo)
    fileName = strcat('./OriginalData/OnlyXsens/WithExo/', dinfo(ind_file).name);
    [exp_frames, discarted_frames] = ObtainXsensData(fileName, 100, false);
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
    
    with_exo_data(fileName) = exp_data;
end

if (transform)
    save(".\DataSets\OnlyXsens\transTrainData", "with_exo_data", "without_exo_data")
else
    save(".\DataSets\OnlyXsens\noTransTrainData", "with_exo_data", "without_exo_data")
end

% save(".\DataSets\OnlyXsens\transTrainData", "with_exo_data", "without_exo_data")
% save(".\DataSets\OnlyXsens\noTransTrainData", "with_exo_data", "without_exo_data")

%% Plot of the remaining variables 

n_exps_no_exo = length(without_exo_data.keys);
n_exps_exo = length(with_exo_data.keys);

n_max_exps = max(n_exps_exo,n_exps_no_exo);
exp_number = 0;
without_exo_keys = without_exo_data.keys;
for i_key = 1:length(without_exo_keys)
    no_exo_key = without_exo_keys(i_key);
    figure_number = 1;
    sub_plot_numb = 1 + (exp_number*2);

    exp_data = without_exo_data(no_exo_key);
    mod_exp_data = exp_data.modified_frames;

    seg_names = fieldnames(mod_exp_data.segments);
    for i_seg=1:length(seg_names)
        sg_name = seg_names{i_seg};
        
        labels = {'x','y','z'};
        labels_quat = {'q1','q2','q3','q4'};

        varnames = fieldnames(mod_exp_data.segments.(sg_name));
        for k_pos=1:length(varnames)
            varname = varnames{k_pos};
            if not(varname=="position")
                continue
            end
            dims_var = size(mod_exp_data.segments.(sg_name).(varname));
            for j_pos=1:dims_var(2)
                figure(figure_number)
                subplot(n_max_exps, 2, sub_plot_numb)
                plot(mod_exp_data.segments.(sg_name).(varname)(:,j_pos))
                figure_number = figure_number + 1;
                if dims_var(2)==3
                    title(strcat(varname, " ", labels{j_pos} ,": ",sg_name))
                else
                    title(strcat(varname, " ", labels_quat{j_pos} ,": ",sg_name))
                end
            end
        end
    end
    jnt_names = fieldnames(mod_exp_data.joints);
    for i_jnt=1:length(jnt_names)
        
        jnt_name = jnt_names{i_jnt};
        labels = {'x','y','z'};
        dims_var = size(mod_exp_data.joints.(jnt_name));
        for j_pos=1:dims_var(2)
            figure(figure_number)
            subplot(n_max_exps, 2, sub_plot_numb)
            plot(mod_exp_data.joints.(jnt_name)(:,j_pos))
            figure_number = figure_number + 1;
            title(strcat(varname, " ", labels{j_pos} ,": ",jnt_name))
        end
    end
    exp_number = exp_number+1;
end


exp_number = 1;
with_exo_keys = with_exo_data.keys;
for i_key = 1:length(with_exo_keys)
    exo_key = with_exo_keys(i_key);
    figure_number = 1;
    sub_plot_numb = (exp_number*2);

    exp_data = with_exo_data(exo_key);
    mod_exp_data = exp_data.modified_frames;

    seg_names = fieldnames(mod_exp_data.segments);
    for i_seg=1:length(seg_names)
        sg_name = seg_names{i_seg};
        
        labels = {'x','y','z'};
        labels_quat = {'q1','q2','q3','q4'};

        varnames = fieldnames(mod_exp_data.segments.(sg_name));
        for k_pos=1:length(varnames)
            varname = varnames{k_pos};
            if not(varname=="position")
                continue
            end
            dims_var = size(mod_exp_data.segments.(sg_name).(varname));
            for j_pos=1:dims_var(2)
                figure(figure_number)
                subplot(n_max_exps, 2, sub_plot_numb)
                plot(mod_exp_data.segments.(sg_name).(varname)(:,j_pos))
                figure_number = figure_number + 1;
                if dims_var(2)==3
                    title(strcat(varname, " ", labels{j_pos} ,": ",sg_name))
                else
                    title(strcat(varname, " ", labels_quat{j_pos} ,": ",sg_name))
                end
            end
        end
    end

    jnt_names = fieldnames(mod_exp_data.joints);
    for i_jnt=1:length(jnt_names)
        
        jnt_name = jnt_names{i_jnt};
        labels = {'x','y','z'};
        dims_var = size(mod_exp_data.joints.(jnt_name));
        for j_pos=1:dims_var(2)
            figure(figure_number)
            subplot(n_max_exps, 2, sub_plot_numb)
            plot(mod_exp_data.joints.(jnt_name)(:,j_pos))
            figure_number = figure_number + 1;
            title(strcat(varname, " ", labels{j_pos} ,": ",jnt_name))
        end
    end
    exp_number = exp_number+1;
    
end
load("DataSets\OnlyXsens\transTrainData.mat");

n_exps_no_exo = length(without_exo_data.keys);
n_exps_exo = length(with_exo_data.keys);
n_max_exps = max(n_exps_exo,n_exps_no_exo);
without_exo_keys = without_exo_data.keys;

fr = 100;
exp_number = 0;
for key_id=1:length(without_exo_keys)
    figure_number = 1;
    sub_plot_numb = 1 + (exp_number*2);

    key = without_exo_keys(key_id);
    DB = without_exo_data(key).unified_db;
    dims = size(DB);
    f = fr*(0:dims(1)/2 )/dims(1);
    for i=1:dims(2)
        data=DB(:,i);
        four = fft(data);
        four = abs(four/length(four));
        four = four(1:(length(four)/2+1));
        four(2:end-1) = 2*four(2:end-1);

        figure(figure_number)
        subplot(n_max_exps, 2, sub_plot_numb)
        plot(f, four)
        figure_number = figure_number + 1;
    end
    exp_number=exp_number+1;
end


with_exo_keys = with_exo_data.keys;
exp_number = 1;

for key_id=1:length(with_exo_keys)
    figure_number = 1;
    sub_plot_numb = (exp_number*2);

    key = with_exo_keys(key_id);
    DB = with_exo_data(key).unified_db;
    dims = size(DB);
    f = fr*(0:dims(1)/2 )/dims(1);
    for i=1:dims(2)
        data=DB(:,i);
        four = fft(data);
        four = abs(four/length(four));
        four = four(1:(length(four)/2+1));
        four(2:end-1) = 2*four(2:end-1);

        figure(figure_number)
        subplot(n_max_exps, 2, sub_plot_numb)
        plot(f, four)
        figure_number = figure_number + 1;
    end
    exp_number=exp_number+1;
end

% nq_freq = fr/2;
% cutoff_low = 0.7;
% 
% [b,a] = butter(4, cutoff_low/nq_freq,"low");
% 


%% 

load("DataSets\OnlyXsens\transTrainData.mat");

n_exps_no_exo = length(without_exo_data.keys);
n_exps_exo = length(with_exo_data.keys);

n_max_exps = max(n_exps_exo,n_exps_no_exo);

without_exo_keys = without_exo_data.keys;

exp_number = 0;

for key_id=1:length(without_exo_keys)
    figure_number = 1;
    sub_plot_numb = 1 + (exp_number*2);

    key = without_exo_keys(key_id);
    DB = without_exo_data(key).unified_db;
    DB_fil = without_exo_data(key).filtered_db;
    
    dims = size(DB);

    for i=1:dims(2)
        data=DB(:,i);
        data_fil = DB_fil(:,i);
        figure(figure_number)
        subplot(n_max_exps, 2, sub_plot_numb)
        plot(data)
        hold on
        plot(data_fil)
        figure_number = figure_number + 1;
    end
    exp_number = exp_number+1;
end

with_exo_keys = with_exo_data.keys;
exp_number = 1;

for key_id=1:length(with_exo_keys)
    figure_number = 1;
    sub_plot_numb = (exp_number*2);
    
    key = with_exo_keys(key_id);
    DB = with_exo_data(key).unified_db;
    DB_fil = with_exo_data(key).filtered_db;

    dims = size(DB);

    for i=1:dims(2)
        data=DB(:,i);
        data_fil = DB_fil(:,i);
        figure(figure_number)
        subplot(n_max_exps, 2, sub_plot_numb)
        plot(data)
        hold on
        plot(data_fil)
        figure_number = figure_number + 1;
    end

    exp_number = exp_number+1;
end
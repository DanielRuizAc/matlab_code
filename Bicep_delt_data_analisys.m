clear

load("DataSets\clustering_data.mat");
load("DataSets\withExoEMG.mat");
load("DataSets\withoutExoEMG.mat");
% load("DataSets\maxEMGVals.mat");

fr = 1000;
nq_freq = fr/2;
cutoff_low = 10;
cutoff_high = 150;

max_emg_vals = struct;

[b,a] = butter(4,[cutoff_low,cutoff_high]/nq_freq,"bandpass");
exp_names = emg_without_exo_data.keys();

channels_muscles=["Right deltoid", "Right Biceps", "Left Deltoid", "Left Biceps"];

for i=1:length(exp_names)
    exp_names{i}
    emg_exp_data = emg_without_exo_data(exp_names{i});
    clust_data = clustering_without_exo(exp_names{i});
    valid_indexes = find((clust_data.ids==2));
    recorted_data = struct;
    channels_data = emg_exp_data.ch_data;
    WinLen = 4000;                                            % Window Length For RMS Calculation
    index=1;
    step=400;
    fieldnames_emg = fieldnames(channels_data);
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        figure(2*i -1)
        subplot(4,1,j)
        plot(ch_data(1,:)*0.001, ch_data(2,:))
        ylim([-1.6e-3, 1.6e-3]);
        xlabel("time (s)")
        ylabel("Sensor signal (V)")
        title(strcat("Sensor output ",channels_muscles(j)))

        f = fr*(0:length(ch_data(2,:))/2)/length(ch_data(2,:));
        four = fft(ch_data(2,:));
        four = abs(four/length(four));
        four = four(1:(length(four)/2+1));
        four(2:end-1) = 2*four(2:end-1);
        figure(2*i)
        subplot(4,1,j)
        plot(f, four)
        xlabel("Frequency Hz")
        ylabel("|fft|")
        title(strcat("Fourier transform ",channels_muscles(j)," data"))
        ch_data(2,:) = filtfilt(b,a, ch_data(2,:));
        if (isfield(max_emg_vals, ch_name)==false)
            max_emg_vals.(ch_name) = 0;
        end
        max_emg_vals.(ch_name) = max(max_emg_vals.(ch_name), max(ch_data(2,:)));
        % ch_data(2,:) = ch_data(2,:)./max_emg_vals.(ch_name);
        emg_indexes = round(ch_data(1,:)/10);
        recorted_data.(ch_name).recorted = ch_data(2, ismember(emg_indexes, valid_indexes));
    end
    emg_without_exo_data(exp_names{i}).recorted = recorted_data;
end


exp_names = emg_with_exo_data.keys();
for i=1:length(exp_names)
    exp_data = emg_with_exo_data(exp_names{i});
    clust_data = clustering_with_exo(exp_names{i});

    valid_indexes = find((clust_data.ids==2));
    recorted_data = struct;
    channels_data = emg_exp_data.ch_data;
    WinLen = 4000;                                            % Window Length For RMS Calculation
    index=1;
    step=400;
    fieldnames_emg = fieldnames(channels_data);
    figure
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        
        figure(2*(i+3) -1)
        subplot(4,1,j)
        plot(ch_data(1,:)*0.001, ch_data(2,:))
        ylim([-1.6e-3, 1.6e-3]);
        xlabel("time (s)")
        ylabel("Sensor signal (V)")
        title(strcat("Sensor output ",channels_muscles(j)))

        f = fr*(0:length(ch_data(2,:))/2)/length(ch_data(2,:));
        four = fft(ch_data(2,:));
        four = abs(four/length(four));
        four = four(1:(length(four)/2+1));
        four(2:end-1) = 2*four(2:end-1);
        figure(2*(i+3))
        subplot(4,1,j)
        plot(f, four)
        xlabel("Frequency Hz")
        ylabel("|fft|")
        title(strcat("Fourier transform ",channels_muscles(j)," data"))
        ch_data(2,:) = filtfilt(b,a, ch_data(2,:));
        if (isfield(max_emg_vals, ch_name)==false)
            max_emg_vals.(ch_name) = 0;
        end
        max_emg_vals.(ch_name) = max(max_emg_vals.(ch_name), max(ch_data(2,:)));
        % ch_data(2,:) = ch_data(2,:)./max_emg_vals.(ch_name);
        emg_indexes = round(ch_data(1,:)/10);
        recorted_data.(ch_name).recorted = ch_data(2, ismember(emg_indexes, valid_indexes));
    end
    emg_with_exo_data(exp_names{i}).recorted = recorted_data;
end


%% 

exp_names = emg_without_exo_data.keys();
for i=1:length(exp_names)
    exp_names{i}
    emg_exp_data = emg_without_exo_data(exp_names{i});
    clust_data = clustering_without_exo(exp_names{i});
    valid_indexes = find((clust_data.ids==2));
    recorted_data = struct;
    channels_data = emg_exp_data.ch_data;
    WinLen = 4000;                                            % Window Length For RMS Calculation
    index=1;
    step=400;
    fieldnames_emg = fieldnames(channels_data);
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        ch_data(2,:) = filtfilt(b,a, ch_data(2,:));
        ch_data(2,:) = ch_data(2,:)./max_emg_vals.(ch_name);
        emg_indexes = round(ch_data(1,:)/10);
        recorted_data.(ch_name).recorted = ch_data(2, ismember(emg_indexes, valid_indexes));
        window_indexes = 1:step:length(recorted_data.(ch_name).recorted)-WinLen;

        rmsv = zeros(1, length(window_indexes));
        arvv = zeros(1, length(window_indexes));
        md_freq = zeros(1, length(window_indexes));
        mean_freq = zeros(1, length(window_indexes));
        index=1;
        for k=window_indexes
            signal_cutted = recorted_data.(ch_name).recorted(k:k+WinLen);
            rmsv(index) = rms(signal_cutted,2);
            arvv(index) = sum(abs(signal_cutted),2)/WinLen;
            md_freq(index) = medfreq(signal_cutted,fr);
            mean_freq(index) = meanfreq(signal_cutted,fr);
            index=index+1;
        end
    
        recorted_data.(ch_name).rmsv = rmsv;
        recorted_data.(ch_name).arvv = arvv;
        recorted_data.(ch_name).md_freq = md_freq;
        recorted_data.(ch_name).mean_freq = mean_freq;

        figure(15)
        subplot(4,1,j)
        hold on
        plot(rmsv)
        ylim([0.05, 0.22])

        figure(16)
        subplot(4,1,j)
        hold on
        plot(arvv)
        ylim([0.05, 0.22])

        figure(17)
        subplot(4,1,j)
        hold on
        plot(md_freq)
        ylim([40, 65])

        figure(18)
        subplot(4,1,j)
        hold on
        plot(mean_freq)
        ylim([40, 65])
    end
    emg_without_exo_data(exp_names{i}).recorted = recorted_data;
end


exp_names = emg_with_exo_data.keys();
for i=1:length(exp_names)
    exp_data = emg_with_exo_data(exp_names{i});
    clust_data = clustering_with_exo(exp_names{i});

    valid_indexes = find((clust_data.ids==2));
    recorted_data = struct;
    channels_data = emg_exp_data.ch_data;
    WinLen = 4000;                                            % Window Length For RMS Calculation
    index=1;
    step=400;
    fieldnames_emg = fieldnames(channels_data);
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        ch_data(2,:) = filtfilt(b,a, ch_data(2,:));
        
        ch_data(2,:) = ch_data(2,:)./max_emg_vals.(ch_name);
        emg_indexes = round(ch_data(1,:)/10);
        recorted_data.(ch_name).recorted = ch_data(2, ismember(emg_indexes, valid_indexes));
        window_indexes = 1:step:length(recorted_data.(ch_name).recorted)-WinLen;

        rmsv = zeros(1, length(window_indexes));
        arvv = zeros(1, length(window_indexes));
        md_freq = zeros(1, length(window_indexes));
        mean_freq = zeros(1, length(window_indexes));
        index=1;
        for k=window_indexes
            signal_cutted = recorted_data.(ch_name).recorted(k:k+WinLen);
            rmsv(index) = rms(signal_cutted,2);
            arvv(index) = sum(abs(signal_cutted),2)/WinLen;
            md_freq(index) = medfreq(signal_cutted,fr);
            mean_freq(index) = meanfreq(signal_cutted,fr);
            index=index+1;
        end
    
        recorted_data.(ch_name).rmsv = rmsv(1:750);
        recorted_data.(ch_name).arvv = arvv(1:750);
        recorted_data.(ch_name).md_freq = md_freq(1:750);
        recorted_data.(ch_name).mean_freq = mean_freq(1:750);
        
        figure(19)
        subplot(4,1,j)
        hold on
        plot(rmsv(1:850))
        ylim([0.05, 0.22])

        figure(20)
        subplot(4,1,j)
        hold on
        plot(arvv(1:850))
        ylim([0.05, 0.22])

        figure(21)
        subplot(4,1,j)
        hold on
        plot(md_freq(1:850))
        ylim([40, 65])

        figure(22)
        subplot(4,1,j)
        hold on
        plot(mean_freq(1:850))
        ylim([40, 65])

    end
    emg_with_exo_data(exp_names{i}).recorted = recorted_data;
end




%%

exp_names = emg_without_exo_data.keys();
for i=1:length(exp_names) 
    exp_names{i}
    emg_exp_data = emg_without_exo_data(exp_names{i});
    channels_data = emg_exp_data.recorted;
    fieldnames_emg = fieldnames(channels_data);
    regression_vals = struct;
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        length_coefs = length(ch_data.rmsv);
        X=[ones(length_coefs,1) (1:length_coefs)'];
        regression_vals.(ch_name).rmsv_reg = X\ch_data.rmsv';
        regression_vals.(ch_name).arvv_reg = X\ch_data.arvv';
        regression_vals.(ch_name).md_freq_reg = X\ch_data.md_freq';
        regression_vals.(ch_name).mean_freq_reg = X\ch_data.mean_freq';
    end
    emg_without_exo_data(exp_names{i}).regression = regression_vals;
end


exp_names = emg_with_exo_data.keys();
for i=1:length(exp_names) 
    exp_names{i}
    emg_exp_data = emg_with_exo_data(exp_names{i});
    channels_data = emg_exp_data.recorted;
    fieldnames_emg = fieldnames(channels_data);
    regression_vals = struct;
    for j=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{j};
        ch_data = channels_data.(ch_name);
        length_coefs = length(ch_data.rmsv);
        X=[ones(length_coefs,1) (1:length_coefs)'];
        regression_vals.(ch_name).rmsv_reg = X\ch_data.rmsv';
        regression_vals.(ch_name).arvv_reg = X\ch_data.arvv';
        regression_vals.(ch_name).md_freq_reg = X\ch_data.md_freq';
        regression_vals.(ch_name).mean_freq_reg = X\ch_data.mean_freq';
    end
    emg_with_exo_data(exp_names{i}).regression = regression_vals;
end



%%

without_exo_experiment= emg_without_exo_data("2023-11-22-131528.txt");
with_exo_experiment= emg_with_exo_data("2023-11-22-135935.txt");

legs = ["Without Exoskeleton", "With Exoskeleton"];

figure
subplot(4,1,1)
plot((1:length(without_exo_experiment.recorted.ch_0.rmsv))*400/1000, ...
    without_exo_experiment.recorted.ch_0.rmsv, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_0.rmsv))*400/1000, ...
    with_exo_experiment.recorted.ch_0.rmsv, "LineWidth", 1.1)
ylim([0.05, 0.25])
xlabel("time (s)")
ylabel("rms %")
legend(legs)
title(strcat("RMS calculated using windows of 4 seconds for ",channels_muscles(1)))


subplot(4,1,2)
plot((1:length(without_exo_experiment.recorted.ch_1.rmsv))*400/1000, ...
    without_exo_experiment.recorted.ch_1.rmsv, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_1.rmsv))*400/1000, ...
    with_exo_experiment.recorted.ch_1.rmsv, "LineWidth", 1.1)
ylim([0.05, 0.25])
xlabel("time (s)")
ylabel("rms %")
legend(legs)
title(strcat("RMS calculated using windows of 4 seconds for ",channels_muscles(2)))


subplot(4,1,3)
plot((1:length(without_exo_experiment.recorted.ch_2.rmsv))*400/1000, ...
    without_exo_experiment.recorted.ch_2.rmsv, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_2.rmsv))*400/1000, ...
    with_exo_experiment.recorted.ch_2.rmsv, "LineWidth", 1.1)
ylim([0.05, 0.25])
xlabel("time (s)")
ylabel("rms %")
legend(legs)
title(strcat("RMS calculated using windows of 4 seconds for ",channels_muscles(3)))

subplot(4,1,4)
plot((1:length(without_exo_experiment.recorted.ch_0.rmsv))*400/1000, ...
    without_exo_experiment.recorted.ch_3.rmsv, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_0.rmsv))*400/1000, ...
    with_exo_experiment.recorted.ch_3.rmsv, "LineWidth", 1.1)
ylim([0.05, 0.25])
xlabel("time (s)")
ylabel("rms %")
legend(legs)
title(strcat("RMS calculated using windows of 4 seconds for ",channels_muscles(4)))

figure
subplot(4,1,1)
plot((1:length(without_exo_experiment.recorted.ch_0.md_freq))*400/1000, ...
    without_exo_experiment.recorted.ch_0.md_freq, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_0.md_freq))*400/1000, ...
    with_exo_experiment.recorted.ch_0.md_freq, "LineWidth", 1.1)
ylim([38, 60])
xlabel("time (s)")
ylabel("Frequency (Hz)")
legend(legs)
title(strcat("Median frequency calculated using windows of 4 seconds for ",channels_muscles(1)))

subplot(4,1,2)
plot((1:length(without_exo_experiment.recorted.ch_1.md_freq))*400/1000, ...
    without_exo_experiment.recorted.ch_1.md_freq, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_1.md_freq))*400/1000, ...
    with_exo_experiment.recorted.ch_1.md_freq, "LineWidth", 1.1)
ylim([36, 55])
xlabel("time (s)")
ylabel("Frequency (Hz)")
legend(legs)
title(strcat("Median frequency calculated using windows of 4 seconds for ",channels_muscles(2)))

subplot(4,1,3)
plot((1:length(without_exo_experiment.recorted.ch_2.md_freq))*400/1000, ...
    without_exo_experiment.recorted.ch_2.md_freq, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_2.md_freq))*400/1000, ...
    with_exo_experiment.recorted.ch_2.md_freq, "LineWidth", 1.1)
ylim([38, 60])
xlabel("time (s)")
ylabel("Frequency (Hz)")
legend(legs)
title(strcat("Median frequency calculated using windows of 4 seconds for ",channels_muscles(3)))

subplot(4,1,4)
plot((1:length(without_exo_experiment.recorted.ch_3.md_freq))*400/1000, ...
    without_exo_experiment.recorted.ch_3.md_freq, "LineWidth", 1.1)
hold on
plot((1:length(with_exo_experiment.recorted.ch_3.md_freq))*400/1000, ...
    with_exo_experiment.recorted.ch_3.md_freq, "LineWidth", 1.1)
ylim([36, 55])
xlabel("time (s)")
ylabel("Frequency (Hz)")
legend(legs)
title(strcat("Median frequency calculated using windows of 4 seconds for ",channels_muscles(3)))
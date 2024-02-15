clear
load("DataSets\OnlyEMG\EMG_DATA.mat")
dinfo = dir('./OriginalData/XSensEmg/WithoutExo/*.txt');
without_exo_data = without_exo_data_emg;% dictionary;

fr = 1000;
nq_freq = fr/2;
cutoff_low = 10;
cutoff_high = 200;

[b,a] = butter(4,[cutoff_low,cutoff_high]/nq_freq,"bandpass");

for ind_file = 1 : length(dinfo)
    if ind_file>1
        break
    end
    emg_exp_data = struct;
    fileName = strcat('./OriginalData/XSensEmg/WithoutExo/', dinfo(ind_file).name)
    % if (without_exo_data.isKey(fileName))
    %     continue
    % end
    [channels_data] = ObtainEMGData(fileName, true);
    fieldnames_emg = fieldnames(channels_data);
    emg_exp_data.ch_data = channels_data;
    for i=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{i};
        channels_data.(ch_name)(2,:) = filtfilt(b,a,abs(channels_data.(ch_name)(2,:)));
    end
    % channels_data(2,:) = abs(filtfilt(b,a,channels_data(2,:)));
    emg_exp_data.rect_filt = channels_data;
    without_exo_data(fileName) = emg_exp_data;
end


dinfo = dir('./OriginalData/XSensEmg/WithExo/*.txt');
with_exo_data = without_exo_data_emg; % dictionary;

fr = 1000;
nq_freq = fr/2;
cutoff_low = 10;
cutoff_high = 200;

[b,a] = butter(4,[cutoff_low,cutoff_high]/nq_freq,"bandpass");

for ind_file = 1 : length(dinfo)
    break;
    emg_exp_data = struct;
    fileName = strcat('./OriginalData/XSensEmg/WithExo/', dinfo(ind_file).name)
    % if (with_exo_data.isKey(fileName))
    %     continue
    % end
    [channels_data] = ObtainEMGData(fileName, true);
    fieldnames_emg = fieldnames(channels_data);
    emg_exp_data.ch_data = channels_data;
    for i=1:length(fieldnames_emg)
        ch_name = fieldnames_emg{i};
        channels_data.(ch_name)(2,:) = filtfilt(b,a,abs(channels_data.(ch_name)(2,:)));
    end
    % channels_data(2,:) = abs(filtfilt(b,a,channels_data(2,:)));
    emg_exp_data.rect_filt = channels_data;
    with_exo_data(fileName) = emg_exp_data;
end

%% 

without_exo_data_emg = without_exo_data;
with_exo_data_emg = with_exo_data;
save("DataSets\OnlyEMG\EMG_DATA", "with_exo_data_emg","without_exo_data_emg")



% %%
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
% 
% 
% with_exo_keys = with_exo_data.keys;
% for key_id=1:length(with_exo_keys)
%     key = with_exo_keys(key_id);
%     filt_ch_data = with_exo_data(key).ch_data; % without_exo_data(key).rect_filt;
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
%     with_exo_data(key).fat_data = fat_data;
% end
% 
% %%
% 
% exp_number = 0;
% without_exo_keys = without_exo_data.keys;
% 
% figure_number = 1;
% for key_id=1:length(without_exo_keys)
%     % if (key_id)>1
%     %     break
%     % end
% 
%     key = without_exo_keys(key_id);
%     fat_data = without_exo_data(key).fat_data;
% 
%     ch_names = fieldnames(fat_data);
% 
%     for ch_i=1:length(ch_names)
%         ch_name = ch_names{ch_i};
%         ch_fat_data = fat_data.(ch_name);
%         disp(ch_name)
%         figure(figure_number + mod(ch_i-1,3))
%         if ch_i<=3
%             subplot(2,3,1)
%             plot(ch_fat_data.rms_dat)
%             title(ch_name, " RMS")
% 
%             subplot(2,3,2)
%             plot(ch_fat_data.md_freq)
%             title(ch_name, " median frequency")
% 
%             subplot(2,3,3)
%             plot(ch_fat_data.mn_freq)
%             title(ch_name, " mean frequency")
%         else
%             subplot(2,3,4)
%             plot(ch_fat_data.rms_dat)
%             title(ch_name, " RMS")
% 
%             subplot(2,3,5)
%             plot(ch_fat_data.md_freq)
%             title(ch_name, " median frequency")
% 
%             subplot(2,3,6)
%             plot(ch_fat_data.mn_freq)
%             title(ch_name, " mean frequency")
%         end
%     end
%     figure_number = figure_number + 3;
% end
% 
% 
% %%
% 
% with_exo_keys = with_exo_data.keys;
% 
% 
% for key_id=1:length(with_exo_keys)
%     % if (key_id)>1
%     %     break
%     % end
% 
%     key = with_exo_keys(key_id);
%     fat_data = with_exo_data(key).fat_data;
% 
%     ch_names = fieldnames(fat_data);
% 
%     for ch_i=1:length(ch_names)
%         ch_name = ch_names{ch_i};
%         ch_fat_data = fat_data.(ch_name);
%         disp(ch_name)
%         figure(figure_number + mod(ch_i-1,3))
%         if ch_i<=3
%             subplot(2,3,1)
%             plot(ch_fat_data.rms_dat)
% 
%             subplot(2,3,2)
%             plot(ch_fat_data.md_freq)
% 
%             subplot(2,3,3)
%             plot(ch_fat_data.mn_freq)
%         else
%             subplot(2,3,4)
%             plot(ch_fat_data.rms_dat)
% 
%             subplot(2,3,5)
%             plot(ch_fat_data.md_freq)
% 
%             subplot(2,3,6)
%             plot(ch_fat_data.mn_freq)
%         end
%     end
%     figure_number = figure_number + 3;
% end
% 
% 
% %% 
% 
% without_exo_data_emg = without_exo_data;
% with_exo_data_emg = with_exo_data;
% save("DataSets\OnlyEMG\EMG_DATA", "with_exo_data_emg","without_exo_data_emg")
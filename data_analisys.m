clear

dinfo = dir('./without exo/*.txt');
without_exo_data = dictionary;

for ind_file = 1 : length(dinfo)
    fileName = strcat('./without exo/', dinfo(ind_file).name);  %just the name
    str_emg = fileread(fileName); % dedicated for reading files as text
    json_x_sens_msgs_emg = splitlines(str_emg); % Splitting the complete text in lines
    
    index_xsens = 1;
    channels_data = {};

    channels_data = {};

    for i=1:(length(json_x_sens_msgs_emg)-1)
        json_packet= jsondecode(json_x_sens_msgs_emg{i});
        bts_data_pkt = json_packet.bts.data;
        frames = json_packet.xsens.data;
        if (isempty(bts_data_pkt)) 
            continue
        end
        for j=1:length(bts_data_pkt)
            ch_data = bts_data_pkt(j).data;
            ch_values = [];
            ch_index = [];
            for k=1:length(ch_data)
                ch_values(k) = ch_data(k).value;
                ch_index(k) = ch_data(k).index;
            end
            if (j>length(channels_data))
                channels_data{j} = [ch_index; ch_values];
            else
                channels_data{j} =[channels_data{j} [ch_index; ch_values]];
            end
        end
        
        for j=1:(length(frames))
            xsens_data(index_xsens) = frames(j);
            index_xsens=index_xsens+1;
        end
    end

    exp_data.bts=channels_data;
    exp_data.xsens=xsens_data;
    without_exo_data(fileName) = exp_data;
end

%%
% fileName = './without exo/2023-10-27-11394.txt';
% str_emg = fileread(fileName); % dedicated for reading files as text
% json_x_sens_msgs_emg = splitlines(str_emg); % Splitting the complete text in lines
% 
% index_bts = 0;
% index_xsens = 1;
% 
% channels_data = {};
% 
% for i=1:(length(json_x_sens_msgs_emg)-1)
%     json_packet= jsondecode(json_x_sens_msgs_emg{i});
%     bts_data_pkt = json_packet.bts.data;
%     frames = json_packet.xsens.data;
%     if (isempty(bts_data_pkt)) 
%         continue
%     end
%     for j=1:length(bts_data_pkt)
%         ch_data = bts_data_pkt(j).data;
%         ch_values = [];
%         ch_index = [];
%         for k=1:length(ch_data)
%             ch_values(k) = ch_data(k).value;
%             ch_index(k) = ch_data(k).index;
%         end
%         if (j>length(channels_data))
%             channels_data{j} = [ch_index; ch_values];
%         else
%             channels_data{j} =[channels_data{j} [ch_index; ch_values]];
%         end
%     end
% 
%     for j=1:(length(frames))
%         bts_data(index_xsens) = frames(j);
%         index_xsens=index_xsens+1;
%     end
% end
% %%
% 
% rightFore = zeros(length(bts_data),1);
% 
% for i=1:length(rightFore)
%     rightFore(i) = bts_data(i).segments.RightForeArm.position.z;
% end
% %% 
% 
% figure
% for i=1:length(channels_data)
%     subplot(length(channels_data),1,i)
%     plot(channels_data{i}(1,:)*0.001,channels_data{i}(2,:))
% end

%%
keys = without_exo_data.keys;
for i=1:length(without_exo_data.keys)
    disp(i)
    key = keys(i);
    exp_data = without_exo_data(key);
    bts_data_exp = exp_data.bts;
    figure
    for j=1:length(bts_data_exp)
        subplot(length(bts_data_exp),1,j)
        plot(bts_data_exp{j}(1,:)*0.001,bts_data_exp{j}(2,:))
    end
end


%% 
dinfo = dir('./with exo/*.txt');
with_exo_data = dictionary;

for ind_file = 1 : length(dinfo)
    fileName = strcat('./with exo/', dinfo(ind_file).name);  %just the name
    str_emg = fileread(fileName); % dedicated for reading files as text
    json_x_sens_msgs_emg = splitlines(str_emg); % Splitting the complete text in lines
    
    index_xsens = 1;
    channels_data = {};

    channels_data = {};

    for i=1:(length(json_x_sens_msgs_emg)-1)
        json_packet= jsondecode(json_x_sens_msgs_emg{i});
        bts_data_pkt = json_packet.bts.data;
        frames = json_packet.xsens.data;
        if (isempty(bts_data_pkt)) 
            continue
        end
        for j=1:length(bts_data_pkt)
            ch_data = bts_data_pkt(j).data;
            ch_values = [];
            ch_index = [];
            for k=1:length(ch_data)
                ch_values(k) = ch_data(k).value;
                ch_index(k) = ch_data(k).index;
            end
            if (j>length(channels_data))
                channels_data{j} = [ch_index; ch_values];
            else
                channels_data{j} =[channels_data{j} [ch_index; ch_values]];
            end
        end
        
        for j=1:(length(frames))
            xsens_data(index_xsens) = frames(j);
            index_xsens=index_xsens+1;
        end
    end

    exp_data.bts=channels_data;
    exp_data.xsens=xsens_data;
    with_exo_data(fileName) = exp_data;
end


%%
keys = with_exo_data.keys;
for i=1:length(with_exo_data.keys)
    disp(i)
    key = keys(i);
    exp_data = with_exo_data(key);
    bts_data_exp = exp_data.bts;
    figure
    for j=1:length(bts_data_exp)
        subplot(length(bts_data_exp),1,j)
        plot(bts_data_exp{j}(1,:)*0.001,bts_data_exp{j}(2,:))
    end
end
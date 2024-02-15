function [channels_data] = ObtainEMGData(filename, combined)

    file_data = fileread(filename);  % Reading files as text
    json_msgs = splitlines(file_data); % Splitting the complete text in lines
    
    channels_data = struct;
    for i=1:(length(json_msgs)-1)
        json_packet= jsondecode(json_msgs{i});
        if combined
            bts_data_pkt = json_packet.bts.data;
        else
            bts_data_pkt = json_packet.data;
        end
        if (isempty(bts_data_pkt)) 
            continue
        end
        for j=1:length(bts_data_pkt)
            ch_channel = bts_data_pkt(j).Channel;
            ch_data = bts_data_pkt(j).data;
            ch_values = [];
            ch_index = [];
            for k=1:length(ch_data)
                ch_values(k) = ch_data(k).value;
                ch_index(k) = ch_data(k).index;
            end
            str_channel = strcat('ch_' ,num2str(ch_channel));
            if isfield(channels_data ,str_channel)
                channels_data.(str_channel) = [channels_data.(str_channel) [ch_index; ch_values]];
            else
                channels_data.(str_channel) = [ch_index; ch_values];
            end
            
        end

    end

end
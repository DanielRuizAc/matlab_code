function [experiment_frames, discarted_values] = ObtainXsensData(filename, review_limit, combined)

    file_data = fileread(filename);  % Reading files as text
    json_msgs = splitlines(file_data); % Splitting the complete text in lines
    
    index = 1;
    for i=1:(length(json_msgs)-1)
        json_packet= jsondecode(json_msgs{i});
        if combined
            frames = json_packet.xsens.data;
        else
            frames = json_packet.data;
        end
        
        for j=1:(length(frames))
            experiment_frames(index) = frames(j);
            index=index+1;
        end
    end

    discarted_values=0;
    for i=1:review_limit
        if experiment_frames(i).frame_number >= experiment_frames(i+1).frame_number
            if experiment_frames(i+1).frame_number == experiment_frames(i+2).frame_number
                discarted_values=i+1;
                break
            end
            discarted_values=i;
            break
        end
    end
    
    experiment_frames = experiment_frames(discarted_values+1:end);
end
function [mod_frames] = CreateDBXSens(frames, desired_segments, desired_joints, transform)
    
    % mod_frames = frames;
    undesired_segs = setxor(fieldnames(frames(1).segments), desired_segments);
    undesired_joints = setxor(fieldnames(frames(1).joints), desired_joints);

    valid_segs = fieldnames(rmfield(frames(1).segments, undesired_segs));
    valid_joints = fieldnames(rmfield(frames(1).joints, undesired_joints));
    
    n_frames = length(frames);
    mod_frames = struct;

    mod_frames.frame_number = zeros(n_frames,1);
    for i_seg=1:length(valid_segs)
        sg_name = valid_segs{i_seg};
        mod_frames.segments.(sg_name).position = zeros(n_frames,3);
        mod_frames.segments.(sg_name).velocity = zeros(n_frames,3);
        mod_frames.segments.(sg_name).acceleration = zeros(n_frames,3);
        mod_frames.segments.(sg_name).ang_vel = zeros(n_frames,3);
        mod_frames.segments.(sg_name).ang_accel = zeros(n_frames,3);
        mod_frames.segments.(sg_name).orientation_rpy = zeros(n_frames,3);

        mod_frames.segments.(sg_name).orientation = zeros(n_frames,4);
    end

    for i_jnt=1:length(valid_joints)
        jnt_name = valid_joints{i_jnt};
        mod_frames.joints.(jnt_name) = zeros(n_frames,3);
    end
    
    for j=1:n_frames
        Pelv_or = frames(j).segments.Pelvis.orientation;
        Pelv_pos = frames(j).segments.Pelvis.position;

        ref_quat = [Pelv_or.q1 Pelv_or.q2 Pelv_or.q3 Pelv_or.q4];
        ref_pos = [Pelv_pos.x Pelv_pos.y Pelv_pos.z]';
        
        if transform
            R = quat2rotm(ref_quat);
        else
            R = eye(3);
        end
        % R = quat2rotm(ref_quat);
        T_fr_global_to_pelv = [R' -(R')*ref_pos; 0 0 0 1];
        
        for i_seg=1:length(valid_segs)
            sg_name = valid_segs{i_seg};
            segment_data = frames(j).segments.(sg_name); % getfield(frames.segments, sg_name);

            % disp(fieldnames(segment_data))
            
            segment_pos = [segment_data.position.x segment_data.position.y ...
                            segment_data.position.z 1];
            Transf_pose = T_fr_global_to_pelv*(segment_pos');
            mod_frames.segments.(sg_name).position(j,:) = Transf_pose(1:3)';
            
            try

                segment_vel = [segment_data.velocity.x segment_data.velocity.y ...
                            segment_data.velocity.z];
                Transf_vel = (R')*(segment_vel');
                segment_acc = [segment_data.acceleration.x segment_data.acceleration.y ...
                                segment_data.acceleration.z];
                Transf_acc = (R')*(segment_acc');
                mod_frames.segments.(sg_name).velocity(j,:) = Transf_vel(1:3)';
                mod_frames.segments.(sg_name).acceleration(j,:) = Transf_acc(1:3)';
            catch
                % disp("There is not velocity and acceleration in file")
            end
            
            

            segment_quat = [segment_data.orientation.q1 segment_data.orientation.q2 ...
                            segment_data.orientation.q3 segment_data.orientation.q4];
            segment_R = quat2rotm(segment_quat);
            transf_R = (R')*segment_R;

            try
                segment_ang_vel = [segment_data.ang_vel.x segment_data.ang_vel.y ...
                            segment_data.ang_vel.z];
                Transf_ang_vel = (R')*(segment_ang_vel');
                segment_ang_acc = [segment_data.ang_accel.x segment_data.ang_accel.y ...
                                segment_data.ang_accel.z];
                Transf_ang_acc = (R')*(segment_ang_acc');
                mod_frames.segments.(sg_name).ang_vel(j,:) = Transf_ang_vel(1:3)';
                mod_frames.segments.(sg_name).ang_accel(j,:) = Transf_ang_acc(1:3)';
            catch
                % disp("There is no angular velocity and acceleration in the file")
            end
            % segment_ang_vel = [segment_data.ang_vel.x segment_data.ang_vel.y ...
            %                 segment_data.ang_vel.z];
            % Transf_ang_vel = (R')*(segment_ang_vel');
            % segment_ang_acc = [segment_data.ang_accel.x segment_data.ang_accel.y ...
            %                 segment_data.ang_accel.z];
            % Transf_ang_acc = (R')*(segment_ang_acc');
            % 
            % mod_frames.segments.(sg_name).position(j,:) = Transf_pose(1:3)';
            % mod_frames.segments.(sg_name).velocity(j,:) = Transf_vel(1:3)';
            % mod_frames.segments.(sg_name).acceleration(j,:) = Transf_acc(1:3)';
            % mod_frames.segments.(sg_name).ang_vel(j,:) = Transf_ang_vel(1:3)';
            % mod_frames.segments.(sg_name).ang_accel(j,:) = Transf_ang_acc(1:3)';

            mod_frames.segments.(sg_name).orientation(j,:) = rotm2quat(transf_R);
            mod_frames.segments.(sg_name).orientation_rpy(j,:) = rotm2eul(transf_R);
        end

        for i_jnt=1:length(valid_joints)
            jnt_name = valid_joints{i_jnt};
            joint_data = frames(j).joints.(jnt_name);
            mod_frames.joints.(jnt_name)(j,:) = [joint_data.x joint_data.y joint_data.z];
        end

    end

    % for j=1:(length(frames))
    % 
    %     Pelv_or = frames(j).segments.Pelvis.orientation;
    %     Pelv_pos = frames(j).segments.Pelvis.position;
    % 
    %     ref_quat = [Pelv_or.q1 Pelv_or.q2 Pelv_or.q3 Pelv_or.q4];
    %     ref_pos = [Pelv_pos.x Pelv_pos.y Pelv_pos.z]';
    % 
    %     R = quat2rotm(ref_quat);
    %     T_fr_global_to_pelv = [R' -(R')*ref_pos; 0 0 0 1];
    % 
    %     if not(isa(class(undesired_segs), 'cell'))
    %         undesired_segs = setxor(fieldnames(frames(j).segments), desired_segments);
    %         undesired_joints = setxor(fieldnames(frames(j).joints), desired_joints);
    %     end 
    % 
    %     mod_frames(j).segments =  rmfield(frames(j).segments, undesired_segs);
    %     mod_frames(j).joints =  rmfield(frames(j).joints, undesired_joints);
    % 
    %     segment_names = fieldnames(mod_frames(j).segments);
    %     joint_names = fieldnames(mod_frames(j).joints);
    % 
    %     n_segments = length(segment_names);
    %     n_joints = length(joint_names);
    %     for i_seg=1:n_segments
    %         sg_name = segment_names{i_seg};
    %         segment_data = frames(j).segments.(sg_name); % getfield(frames.segments, sg_name);
    % 
    %         segment_pos = [segment_data.position.x segment_data.position.y ...
    %                         segment_data.position.z 1];
    %         Transf_pose = T_fr_global_to_pelv*(segment_pos');
    %         segment_vel = [segment_data.velocity.x segment_data.velocity.y ...
    %                         segment_data.velocity.z];
    %         Transf_vel = (R')*(segment_vel');
    %         segment_acc = [segment_data.acceleration.x segment_data.acceleration.y ...
    %                         segment_data.acceleration.z];
    %         Transf_acc = (R')*(segment_acc');
    % 
    %         segment_quat = [segment_data.orientation.q1 segment_data.orientation.q2 ...
    %                         segment_data.orientation.q3 segment_data.orientation.q4];
    %         segment_R = quat2rotm(segment_quat);
    %         transf_R = (R')*segment_R;
    % 
    %         segment_ang_vel = [segment_data.ang_vel.x segment_data.ang_vel.y ...
    %                         segment_data.ang_vel.z];
    %         Transf_ang_vel = (R')*(segment_ang_vel');
    %         segment_ang_acc = [segment_data.ang_accel.x segment_data.ang_accel.y ...
    %                         segment_data.ang_accel.z];
    %         Transf_ang_acc = (R')*(segment_ang_acc');
    % 
    %         mod_frames(j).segments.(sg_name).position = Transf_pose(1:3)';
    %         mod_frames(j).segments.(sg_name).velocity = Transf_vel(1:3)';
    %         mod_frames(j).segments.(sg_name).acceleration = Transf_acc(1:3)';
    %         mod_frames(j).segments.(sg_name).ang_vel = Transf_ang_vel(1:3)';
    %         mod_frames(j).segments.(sg_name).ang_accel = Transf_ang_acc(1:3)';
    % 
    %         mod_frames(j).segments.(sg_name).orientation = rotm2quat(transf_R);
    %         mod_frames(j).segments.(sg_name).orientation_rpy = rotm2eul(transf_R);
    %     end
    %     for i_jnt=1:n_joints
    %         jnt_name = joint_names{i_jnt};
    %         joint_data = frames(j).joints.(jnt_name); % getfield(frames.joints, jnt_name);
    % 
    %         mod_frames(j).joints.(jnt_name) = [joint_data.x joint_data.y joint_data.z]; % setfield(mod_frames.joints, jnt_name, [joint_data.x joint_data.y joint_data.z]);
    %     end
    % end
    % 
end


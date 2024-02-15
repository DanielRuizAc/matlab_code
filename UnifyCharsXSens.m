function [unifiedDB] = UnifyCharsXSens(modFrames, desiredValues, desiredDimsJoint)
    unifiedDB = [];
    seg_names = fieldnames(modFrames.segments);
    for i_seg=1:length(seg_names)
        sg_name = seg_names{i_seg};
        varnames = fieldnames(modFrames.segments.(sg_name));
        desiredvars = intersect(varnames, fieldnames(desiredValues));
    
        for k_pos=1:length(desiredvars)
            varname = varnames{k_pos};
            if size(unifiedDB)==zeros()
                unifiedDB = modFrames.segments.(sg_name).(varname)(:, desiredValues.(varname));
            else
                unifiedDB = [unifiedDB modFrames.segments.(sg_name).(varname)(:, desiredValues.(varname))];
            end
        end
    
    
    end
    
    jnt_names = fieldnames(modFrames.joints);
    for i_jnt=1:length(jnt_names)
        jnt_name = jnt_names{i_jnt};
        unifiedDB = [unifiedDB modFrames.joints.(jnt_name)(:, desiredDimsJoint)];
    end
end


function [MonoConnectionsTable, MonoConnCCGs, t] = FindMonosynapticConnections(all_data, minNspikes, nreps, smoothingwin)
%FINDMONOSYNAPTICCONNECTIONS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    smoothingwin = [];
end

groups_vec = {};
recordings_vec = {};
senders_vec = {};
receivers_vec = {};
dcw_vec = [];
EorI_vec = {};
cellTypePairs_vec = {};
layerPairs_vec = {};

MonoConnCCGs = []; % will be [time lags x unit pairs], index columns in parallel with height of MonoConnectionsTable

groupNames = fieldnames(all_data);
for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};
        fprintf('Working on recording %s\n',recName);
        cellIDs = fieldnames(all_data.(groupName).(recName));

        numUnits = length(cellIDs);

        figure('Name',recName);
        tiledlayout('flow');

        for cellID_num_A = 1:numUnits
            cellID_A = cellIDs{cellID_num_A};

            % only need to get the recording duration once
            if cellID_num_A == 1
                duration_ms = ceil(all_data.(groupName).(recName).(cellID_A).Recording_Duration * 1000);
            end

            %IsSingleUnit_A = all_data.(groupName).(recName).(cellID_A).IsSingleUnit;
            ISI_violations_percent_A = all_data.(groupName).(recName).(cellID_A).ISI_violations_percent;
            spikeTimes_A = round(all_data.(groupName).(recName).(cellID_A).SpikeTimes_all / 30);
            
            %if IsSingleUnit_A
            if (ISI_violations_percent_A <= 1) && (length(spikeTimes_A) >= minNspikes)
                cell_type_A = all_data.(groupName).(recName).(cellID_A).Cell_Type;
                channel_A = all_data.(groupName).(recName).(cellID_A).Template_Channel;
                layer_A = all_data.(groupName).(recName).(cellID_A).LaminarLabel;
                for cellID_num_B = 1:length(cellIDs)
                    cellID_B = cellIDs{cellID_num_B};
                    %IsSingleUnit_B = all_data.(groupName).(recName).(cellID_B).IsSingleUnit;
                    ISI_violations_percent_B = all_data.(groupName).(recName).(cellID_B).ISI_violations_percent;
                    channel_B = all_data.(groupName).(recName).(cellID_B).Template_Channel;
                    layer_B = all_data.(groupName).(recName).(cellID_B).LaminarLabel;
                    cell_type_B = all_data.(groupName).(recName).(cellID_B).Cell_Type;
                    spikeTimes_B = round(all_data.(groupName).(recName).(cellID_B).SpikeTimes_all / 30);
                    
                    if ~strcmp(cellID_A,cellID_B) && (ISI_violations_percent_B <= 1) && (length(spikeTimes_B) >= minNspikes) && ~((channel_A==channel_B) && strcmp(cell_type_A,cell_type_B)) % only do CCG if A and B are different units

                        [ccg_i,AB_sim,t] = JitterCorrectedCCG(spikeTimes_A, spikeTimes_B, duration_ms, 50, nreps, smoothingwin);

                        %% Classify significant A->B connection if it exists
                        [~,I_max] = max(abs(ccg_i(t>=-10 & t<=10))); % I_max = 11 is the same as t = 0
                        t_max = I_max - 11;
                        if (t_max > 1) && (ccg_i(t==t_max) > AB_sim(t==t_max,2))
                            sig = 1;
                        elseif (t_max > 1) && (ccg_i(t==t_max) < AB_sim(t==t_max,1))
                            sig = -1;
                        else
                            sig = 0;
                        end

                        if sig ~= 0
                            % plot the significant CCG
                            nexttile;
                            plot(t,AB_sim,"Color","#D95319"); hold on; plot(t,ccg_i,"Color","#0072BD"); xline([-10 10]); hold off;
                            title(gca,strcat(cellID_A,' -> ',cellID_B,' | ',cell_type_A,' -> ',cell_type_B,' | ',layer_A,'->',layer_B));

                            % calculate the directed connection weight (DCW) as in Jia et al., 2021
                            dcw = sum(ccg_i(t>0 & t<13)) - sum(ccg_i(t>-13 & t<0)); % positive if A->B excitatory, negative if B->A excitatory, vice versa for inhibitory

                            % add to table of significant connections
                            groups_vec{end+1,1} = groupName;
                            recordings_vec{end+1,1} = recName;
                            senders_vec{end+1,1} = cellID_A;
                            receivers_vec{end+1,1} = cellID_B;
                            cellTypePairs_vec{end+1,1} = strcat(cell_type_A,'->',cell_type_B);
                            layerPairs_vec{end+1,1} = strcat(layer_A,'->',layer_B);
                            if sig > 0
                                EorI_vec{end+1,1} = 'E';
                            else
                                EorI_vec{end+1,1} = 'I';
                            end
                            dcw_vec(end+1,1) = dcw;

                            MonoConnCCGs(:,end+1) = ccg_i;
                        end

                    end
                end
            end
        end
    end
end

varNames = {'group','recording','sender','receiver','CellTypes','layers','EorI','DCW'};

MonoConnectionsTable = table(groups_vec, recordings_vec, senders_vec, ...
    receivers_vec, cellTypePairs_vec, layerPairs_vec, EorI_vec, dcw_vec, ...
    'VariableNames',varNames);

end
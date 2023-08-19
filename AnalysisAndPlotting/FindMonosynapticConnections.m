function [MonoConnectionsTable, MonoConnCCGs, SimultaneousABs, PointwiseABs, t] = FindMonosynapticConnections(all_data, minNspikes, nreps, smoothingwin)
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
sig_vec = [];
EorI_vec = {};
ExcessSync_vec = [];
cellTypePairs_vec = {};
sender_responsivity_vec = {};
receiver_responsivity_vec = {};
nSpkPairs_vec = [];
layerPairs_vec = {};

MonoConnCCGs = []; % will be [time lags x unit pairs], index columns in parallel with height of MonoConnectionsTable
%SimultaneousABs = []; % will be [time lags x unit pairs x 2], 3rd dim 1 is lower, 2 is upper
SimultaneousABs_lower = [];
SimultaneousABs_upper = [];
%PointwiseABs = [];
PointwiseABs_upper = [];
PointwiseABs_lower = [];

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
            if (ISI_violations_percent_A <= 1.5) && (length(spikeTimes_A) >= minNspikes)
                cell_type_A = all_data.(groupName).(recName).(cellID_A).Cell_Type;
                channel_A = all_data.(groupName).(recName).(cellID_A).Template_Channel;
                layer_A = all_data.(groupName).(recName).(cellID_A).LaminarLabel;
                
                stimResponsivity_A = all_data.(groupName).(recName).(cellID_A).StimResponsivity;
                if stimResponsivity_A == 1
                    stimResponsivity_A = '+';
                elseif stimResponsivity_A == 0
                    stimResponsivity_A = 'nr';
                else
                    stimResponsivity_A = '-';
                end
                
                
                for cellID_num_B = 1:length(cellIDs)
                    cellID_B = cellIDs{cellID_num_B};
                    %IsSingleUnit_B = all_data.(groupName).(recName).(cellID_B).IsSingleUnit;
                    ISI_violations_percent_B = all_data.(groupName).(recName).(cellID_B).ISI_violations_percent;
                    channel_B = all_data.(groupName).(recName).(cellID_B).Template_Channel;
                    layer_B = all_data.(groupName).(recName).(cellID_B).LaminarLabel;
                    cell_type_B = all_data.(groupName).(recName).(cellID_B).Cell_Type;

                    stimResponsivity_B = all_data.(groupName).(recName).(cellID_B).StimResponsivity;
                    if stimResponsivity_B == 1
                        stimResponsivity_B = '+';
                    elseif stimResponsivity_B == 0
                        stimResponsivity_B = 'nr';
                    else
                        stimResponsivity_B = '-';
                    end

                    spikeTimes_B = round(all_data.(groupName).(recName).(cellID_B).SpikeTimes_all / 30);
                    
                    if ~strcmp(cellID_A,cellID_B) && (ISI_violations_percent_B <= 1.5) && (length(spikeTimes_B) >= minNspikes) && ~((channel_A==channel_B) && strcmp(cell_type_A,cell_type_B)) % only do CCG if A and B are different units

                        [ccg_i,AB_sim,AB_pt,t,nSpkPairs] = JitterCorrectedCCG(spikeTimes_A, spikeTimes_B, duration_ms, 50, nreps, smoothingwin);
                        
                        nSpkPairs_vec(end+1,1) = nSpkPairs;

                        %% logical classification schemes
                        % first, check for significant synchrony using
                        % two-sided peak
                        [~,I_M] = max(abs(ccg_i(t>=-12 & t<=12))); % I_M = 13 is where t = 0
                        t_M = I_M - 13;
                        if (t_M==0) && (ccg_i(t==t_M) > AB_pt(t==t_M,2)) % significantly synchronous
                            sig = 1;
                            EorI = 'S';
                        else
                            % use one-side peak to assess monosynaptic
                            % connectivity so that if A and B are
                            % bidirectionally connected and the B->A peak
                            % is stronger, A->B will still be significant
                            [~,I_m] = max(abs(ccg_i(t>=0 & t<=12))); % I_m = 1 is where t = 0
                            t_m = I_m - 1;
                            if (t_m >= 1) && (t_m <= 5) && (ccg_i(t==t_m) > AB_sim(t==t_m,2)) % excitatory monosynaptic
                                sig = 1;
                                EorI = 'E';
                            elseif (t_m >= 1) && (t_m <= 5) && (ccg_i(t==t_m) < AB_sim(t==t_m,1)) % inhibitory monosynaptic
                                sig = 1;
                                EorI = 'I';
                            else % not significant ("NS")
                                sig = 0;
                                EorI = 'NS';
                            end
                        end
                        
                        sig_vec(end+1,1) = sig;
                        EorI_vec{end+1,1} = EorI;

                        if sig ==1
                            % plot the significant CCG
                            nexttile;
                            plot(t,AB_sim,"Color","#D95319"); hold on; plot(t,ccg_i,"Color","#0072BD"); xline([-5 5]); hold off;
                            title(gca,strcat(cellID_A,' -> ',cellID_B,' | ',cell_type_A,' -> ',cell_type_B,' | ',layer_A,'->',layer_B));
                            % title(gca,strcat(EorI,' | ',cellID_A,'->',cellID_B,' | ',cell_type_A,'->',cell_type_B));
                        end

                        % calculate the directed connection weight (DCW) as in Jia et al., 2021
                        dcw_vec(end+1,1) = sum(ccg_i(t>0 & t<13)) - sum(ccg_i(t>-13 & t<0)); % positive if A->B excitatory, negative if B->A excitatory, vice versa for inhibitory?

                        % get excess synchrony
                        ExcessSync_vec(end+1,1) = ccg_i(t==0);

                        % add to table of CCG metrics
                        groups_vec{end+1,1} = groupName;
                        recordings_vec{end+1,1} = recName;
                        senders_vec{end+1,1} = cellID_A;
                        receivers_vec{end+1,1} = cellID_B;
                        cellTypePairs_vec{end+1,1} = strcat(cell_type_A,'->',cell_type_B);

                        sender_responsivity_vec{end+1,1} = stimResponsivity_A;
                        receiver_responsivity_vec{end+1,1} = stimResponsivity_B;

                        layerPairs_vec{end+1,1} = strcat(layer_A,'->',layer_B);

                        MonoConnCCGs(:,end+1) = ccg_i;
                        SimultaneousABs_lower(:,end+1) = AB_sim(:,1);
                        SimultaneousABs_upper(:,end+1) = AB_sim(:,2);
                        PointwiseABs_lower(:,end+1) = AB_pt(:,1);
                        PointwiseABs_upper(:,end+1) = AB_pt(:,2);

                    end
                end
            end
        end
    end
end

PointwiseABs = cat(3, PointwiseABs_lower, PointwiseABs_upper);
SimultaneousABs = cat(3, SimultaneousABs_lower, SimultaneousABs_upper);

% varNames = {'group','recording','sender','receiver','CellTypes','layers','EorI','DCW'};
varNames = {'group','recording','UnitA','UnitB','CellTypes','layers','StimResp_A','StimResp_B','N_SpikePairs','Significance','EorI','DCW','ExcessSynchrony'};

% MonoConnectionsTable = table(groups_vec, recordings_vec, senders_vec, ...
%     receivers_vec, cellTypePairs_vec, layerPairs_vec, EorI_vec, dcw_vec, ...
%     'VariableNames',varNames);

MonoConnectionsTable = table(groups_vec, recordings_vec, senders_vec, ...
    receivers_vec, cellTypePairs_vec, layerPairs_vec, sender_responsivity_vec, ...
    receiver_responsivity_vec, nSpkPairs_vec, sig_vec, EorI_vec, ...
    dcw_vec, ExcessSync_vec, ...
    'VariableNames',varNames);

end
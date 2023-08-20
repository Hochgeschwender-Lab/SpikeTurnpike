function Plot_PSTHs_V2(all_data, relative_time_ms, trialTagsLabels, PSTH_type, excludeMUA, firstOrLast,nTrials)
% This one computes PSTHs on the fly so that the number of trials can be
% specified.
%
% PSTH_type: 'raw', 'conv'
% excludeMUA: 0 to keep units with > 1% ISI violations, or 1 to exclude
%             them
% firstOrLast: 'first' to plot first nTrials trials, 'last' to plot last
%              nTrials trials
% nTrials: integer number of trials to use

%% Construct Gaussian filter for PSTH smoothing
gausssigma=3;
gausswindow=21;
tempfil = exp(- (-(gausswindow-1)/2:(gausswindow-1)/2)'.^2 ./ (2*gausssigma^2));
tempfilter = tempfil/sum(tempfil); % sum(tempfil) == sqrt(2*pi)*gausssigma, so it's the correct denominator

%% Initialize arrays for Gramm input
% all of these should have n_rows = n_units
all_PSTHs = []; % [units x time] matrix of PSTHs
groupsVec = {};
trialTagsVec = {};
cellTypesVec = {}; % elements are e.g. 'FS', 'RS', ...
cellTypesVec2 = {}; % elements are e.g. '+FS', '-FS', '+RS', 'nsRS', ...
responsivityVec = {};
layerVec = {};

groupNames = fieldnames(all_data);
for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};
    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};
        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            %IsSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(mouseName).(cellID).ISI_violations_percent;
            
            %if (excludeMUA && IsSingleUnit) || ~excludeMUA
            if (excludeMUA && (ISI_violations_percent <= 1)) || ~excludeMUA
                for trialTagInd = 1:length(trialTagsLabels)

                    %% Compute PSTHs
                    SpikeTrains_ms = all_data.(groupName).(mouseName).(cellID).SpikeTrains_for_PSTHs{1,trialTagInd};
                    % Narrow down to selected number of trials
                    if nargin > 5
                        if strcmp(firstOrLast,'first')
                            SpikeTrains_ms = SpikeTrains_ms(1:nTrials,:);
                        else
                            SpikeTrains_ms = SpikeTrains_ms(end-nTrials:end,:);
                        end
                    end
                    if strcmp(PSTH_type, 'raw')
                        this_PSTH = mean(SpikeTrains_ms,1) * 1000; % *1000 to convert to Hz (spikes per second)
                    else %strcmp(PSTH_type, 'conv')
                        PSTH_raw = mean(SpikeTrains_ms,1) * 1000; % *1000 to convert to Hz (spikes per second)
                        this_PSTH = conv(PSTH_raw, tempfilter, 'same');
                    end
                    all_PSTHs(end+1,:) = this_PSTH;

                    groupsVec{end+1,1} = groupName;
    
                    cell_type = all_data.(groupName).(mouseName).(cellID).Cell_Type;
                    cellTypesVec{end+1,1} = cell_type;
                    trialTagsVec{end+1,1} = trialTagsLabels{trialTagInd};
%                    layerVec{end+1,1} = all_data.(groupName).(mouseName).(cellID).LaminarLabel;
        
                    responsivityNum = all_data.(groupName).(mouseName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        responsivityVec{end+1,1} = responsivityLabel;
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        responsivityVec{end+1,1} = responsivityLabel;
                    else
                        responsivityLabel = '-';
                        responsivityVec{end+1,1} = responsivityLabel;
                    end
    
                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,cell_type);
                end
            end

        end
    end
end

%% Plotting
figure;
g = gramm('x',relative_time_ms, 'y',all_PSTHs, 'color',groupsVec);
g.facet_grid(cellTypesVec, responsivityVec, "scale","free_y");
%g.facet_grid(layerVec, cellTypesVec, "scale","free_y");
%g.facet_grid([], cellTypesVec, "scale","free_y");
g.fig(trialTagsVec);
g.stat_summary('type','sem', 'setylim','true');
%g.geom_line();
g.set_names('x','Time (ms)', 'y','Instantaneous FR (Hz)', 'Color','', 'Row','', 'Column','');
%g.set_order_options('row',{'SG','L4','IG'});
g.draw;

%% TO DO: Print sample sizes as a table (cell type x modulation)

end


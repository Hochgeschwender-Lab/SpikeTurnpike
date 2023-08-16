function ISIs_CV_against_Mode(all_data, cell_types, trialTagsLabels)
% Plot and compare ISI coefficient of variation (CV) between groups.
%
% INPUTS
% - all_data
% - cell_types: cell array of cell types to plot. e.g., {'MSN','TAN'}
% - plot_points: 0 to only plot bars (+/- SEM), 1 to plot points on top.
%       Useful to show the distribution.

groupNames = fieldnames(all_data);

groupsVec = {};
cellTypesVec = {};
cellTypesVec2 = {};
trialTypesVec = {};
CV_vec = [];
CV_baseline_vec = [];
CV2_vec = [];
ISI_PDF_peak_vec = [];
meanISI_vec = [];
AUROC_vec = [];

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        cellIDs = fieldnames(all_data.(groupName).(recName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            %isSingleUnit = all_data.(groupName).(recName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(recName).(cellID).ISI_violations_percent;
            if any(strcmp(cell_types, thisCellType)) % && (ISI_violations_percent <= 1.5)

                trialTagLabels_numeric = unique(all_data.(groupName).(recName).(cellID).Stim_Intensity);

                for trialTagInd = 1:length(trialTagsLabels)
                    groupsVec{end+1,1} = groupName;
                    cellTypesVec{end+1,1} = thisCellType;
                    trialTypesVec{end+1,1} = trialTagsLabels{trialTagInd};

                    responsivityNum = all_data.(groupName).(recName).(cellID).StimResponsivity;
                    if responsivityNum == 1
                        responsivityLabel = '+';
                        % responsivityVec{end+1,1} = responsivityLabel;
                    elseif responsivityNum == 0
                        responsivityLabel = 'nr';
                        % responsivityVec{end+1,1} = responsivityLabel;
                    else
                        responsivityLabel = '-';
                        % responsivityVec{end+1,1} = responsivityLabel;
                    end
                    cellTypesVec2{end+1,1} = strcat(responsivityLabel,thisCellType);

                    trialInds = find(all_data.(groupName).(recName).(cellID).Stim_Intensity == trialTagLabels_numeric(trialTagInd));
    
                    % Calculate CV of ISIs
                    ISI_vec = cell2mat(cellfun(@diff, all_data.(groupName).(recName).(cellID).SpikeTimes_stim(1,trialInds),'UniformOutput',false)) / 30;
                    % ISI_skew = skewness(ISI_vec);
    
                    CV_vec(end+1,1) = std(ISI_vec) / mean(ISI_vec);

                    CV_baseline_vec(end+1,1) = all_data.(groupName).(recName).(cellID).ISI_baseline_CV;

                    meanISI_vec(end+1,1) = mean(ISI_vec);
    
                    % Calculate CV2 of ISIs, which is more robust to FR
                    % fluctuations
                    [CV2, ~] = CV2ISI_ISI(ISI_vec);
                    CV2_vec(end+1,1) = CV2;
    
                    all_data.(groupName).(recName).(cellID).ISI_CV2(trialTagInd) = CV2;

                    ISI_PDF_peak_vec(end+1,1) = all_data.(groupName).(recName).(cellID).ISI_pdf_peak_xy(1);

                    AUROC_vec(end+1,1) = all_data.(groupName).(recName).(cellID).StimProb(1);
                end
            end
        end
    end
end

%% Plotting
% Fig 1: histogram of ISI PDF peaks (modes)
figure;
g = gramm('x',ISI_PDF_peak_vec, 'color',groupsVec, 'subset',strcmp(trialTypesVec,'Zero'));
g.stat_bin("geom","overlaid_bar", 'normalization','count', 'edges',0:1:max(ISI_PDF_peak_vec));
g.set_names('x','ISI Peak (ms)', "color","");
g.draw;

% Fig 2: CV against ISI peak (similar to Shin & Moore figure 2E)
figure;
g = gramm('x',ISI_PDF_peak_vec, 'y',CV_baseline_vec, 'z',AUROC_vec, 'color',groupsVec, 'subset',strcmp(trialTypesVec,'Zero'));
g.geom_point;
g.set_names('x','ISI Peak (ms)', 'y','CV', 'z','AUROC', 'color','');
g.draw;

figure;
g = gramm('x',CV2_vec, 'y',AUROC_vec, 'color',groupsVec, 'subset',strcmp(trialTypesVec,'Zero'));
g.geom_point;
g.set_names('x','CV2', 'y','AUROC', 'color','');
g.draw;

% Fig 3: histogram of mean ISIs
figure;
g = gramm('x',meanISI_vec, 'color',groupsVec, 'subset',strcmp(trialTypesVec,'Zero'));
g.stat_bin("geom","overlaid_bar", 'normalization','count', 'edges',0:2:max(meanISI_vec));
g.set_names('x','Mean ISI (ms)', "color","");
g.draw;

%% make stats table
if nargout > 0
    data_table_CVs = table(groupsVec, cellTypesVec, trialTypesVec, CV_vec, CV2_vec, 'VariableNames',{'Group','CellType','Condition','CV','CV2'});
end

end


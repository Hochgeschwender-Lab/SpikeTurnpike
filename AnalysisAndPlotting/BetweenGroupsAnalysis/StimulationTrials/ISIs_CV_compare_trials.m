function [data_table_CVs, all_data] = ISIs_CV_compare_trials(all_data, cell_types, trialTagsLabels, plot_points)
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
CV2_vec = [];

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
            if any(strcmp(cell_types, thisCellType)) && (ISI_violations_percent < 1.5)

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
    
                    % Calculate CV2 of ISIs, which is more robust to FR
                    % fluctuations
                    [CV2, ~] = CV2ISI_ISI(ISI_vec);
                    CV2_vec(end+1,1) = CV2;
    
                    all_data.(groupName).(recName).(cellID).ISI_CV2(trialTagInd) = CV2;
                end
            end
        end
    end
end

%% Plotting
% Fig 1: bar plot of ISI CV
figure;
g = gramm('x',groupsVec, 'y',CV_vec, 'color',groupsVec);
g.facet_grid(cellTypesVec2, trialTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','ISI CV', 'Color','', 'Row','', 'Column','');
g.no_legend;
g.draw();
if plot_points
    g.update('x',groupsVec, 'y',CV_vec, "color",groupsVec);
    g.geom_jitter();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
end

figure; % Probability density
g = gramm('x',CV_vec(~isnan(CV_vec)), 'color',groupsVec(~isnan(CV_vec)));
g.facet_grid(cellTypesVec2(~isnan(CV_vec)), trialTypesVec(~isnan(CV_vec)), "scale","independent");
g.stat_density('npoints',1000);
g.set_names('x','ISI CV', 'y','Probability Density', 'color','', 'Row','', 'column','');
g.draw;

% Fig 2: bar plot of CV2
figure;
g = gramm('x',groupsVec, 'y',CV2_vec, 'color',groupsVec);
g.facet_grid(cellTypesVec2, trialTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
%g.stat_boxplot();
g.set_names('x','', 'y','ISI CV2', 'Color','', 'Row','', 'Column','');
g.no_legend;
g.draw();

if plot_points
    g.update('x',groupsVec, 'y',CV2_vec, "color",groupsVec);
    g.geom_jitter();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
end

figure; % Probability density
g = gramm('x',CV2_vec(~isnan(CV2_vec)), 'color',groupsVec(~isnan(CV2_vec)));
g.facet_grid(cellTypesVec2(~isnan(CV2_vec)), trialTypesVec(~isnan(CV2_vec)), "scale","independent");
g.stat_density('npoints',1000);
g.set_names('x','ISI CV2', 'y','Probability Density', 'color','', 'Row','', 'column','');
g.draw;

%% make stats table
if nargout > 0
    data_table_CVs = table(groupsVec, cellTypesVec, trialTypesVec, CV_vec, CV2_vec, 'VariableNames',{'Group','CellType','Condition','CV','CV2'});
end

end


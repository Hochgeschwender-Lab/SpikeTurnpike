function [data_table_CVs, all_data] = ISIs_CV_compare(all_data, cell_types, plot_points)
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
            isSingleUnit = all_data.(groupName).(recName).(cellID).IsSingleUnit;
            if any(strcmp(cell_types, thisCellType)) && isSingleUnit
                groupsVec{end+1,1} = groupName;
                cellTypesVec{end+1,1} = thisCellType;
                % CV_vec(end+1,1) = all_data.(groupName).(mouseName).(cellID).ISI_baseline_CV;

                % Calculate CV of ISIs
                ISI_vec = all_data.(groupName).(recName).(cellID).ISI_baseline_vec;
                % ISI_skew = skewness(ISI_vec);
                % ISI_vec = ISI_vec(ISI_vec <= 0.4*1000);

                CV_vec(end+1,1) = std(ISI_vec) / mean(ISI_vec);

                % Calculate CV2 of ISIs, which is more robust to FR
                % fluctuations
                spikeTimes_samples = all_data.(groupName).(recName).(cellID).SpikeTimes_all;
                sr = all_data.(groupName).(recName).(cellID).Sampling_Frequency;
                spikeTimes_s = spikeTimes_samples/sr;
                
                [CV2, ~] = CV2ISI(spikeTimes_s);
                CV2_vec(end+1,1) = CV2;

                all_data.(groupName).(recName).(cellID).ISI_CV2 = CV2;
            end
        end
    end
end

%% Plotting
% Fig 1: bar plot of ISI CV
figure;
g = gramm('x',groupsVec, 'y',CV_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
g.set_names('x','', 'y','ISI CV', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV', 'png');
if plot_points
    g.update('x',groupsVec, 'y',CV_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',CV_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','ISI CV', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV_PDF', 'png');

% Fig 2: bar plot of CV2
figure;
g = gramm('x',groupsVec, 'y',CV2_vec, 'color',groupsVec);
%g.facet_grid([],cellTypesVec, "scale","independent");
g.stat_summary('type','sem', 'geom',{'bar','black_errorbar'}, 'width',1.2, 'setylim',true);
%g.stat_boxplot();
g.set_names('x','', 'y','ISI CV2', 'Color','', 'Column','');
g.no_legend;
g.draw();
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2', 'png');
if plot_points
    g.update('x',groupsVec, 'y',CV2_vec, "color",groupsVec);
    g.geom_point();
    g.set_color_options('lightness',40);
    g.set_point_options("markers","^", "base_size",3);
    g.no_legend;
    g.draw;
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_withPoints', 'fig');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_withPoints', 'svg');
    saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_withPoints', 'png');
end
figure; % Probability density
g = gramm('x',CV2_vec, 'color',groupsVec);
g.stat_density('npoints',1000);
g.set_names('x','ISI CV2', 'y','Probability Density', 'color','');
g.draw;
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_PDF', 'fig');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_PDF', 'svg');
saveas(gcf, '/home/cresp1el-local/Documents/MATLAB/hd_project_sinda/hd_project_sinda/Figs_and_stats/ISI_CV/ISI_CV2_PDF', 'png');

%% make stats table
data_table_CVs = table(groupsVec, cellTypesVec, CV_vec, CV2_vec, 'VariableNames',{'Group','CellType','CV','CV2'});

end


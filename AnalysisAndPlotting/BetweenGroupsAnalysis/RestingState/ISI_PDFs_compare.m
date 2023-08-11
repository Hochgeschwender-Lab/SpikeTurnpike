function ISI_PDFs_compare(all_data, cell_types)
%ISI_PDFS_COMPARE Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);

groupsVec = {};
cellTypesVec = {};
ISIs = {}; % [units x times]
x = [];
ISI_PDFs_X = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    recNames = fieldnames(all_data.(groupName));

    for recNum = 1:length(recNames)
        recName = recNames{recNum};

        cellIDs = fieldnames(all_data.(groupName).(recName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(recName).(cellID).Cell_Type;
            %isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            ISI_violations_percent = all_data.(groupName).(recName).(cellID).ISI_violations_percent;
            if any(strcmp(thisCellType, cell_types)) && (ISI_violations_percent < 1)
                groupsVec{end+1} = groupName;
                cellTypesVec{end+1,1} = thisCellType;
                ISI_vec = all_data.(groupName).(recName).(cellID).ISI_baseline_vec;
                ISIs{end+1} = log10(ISI_vec);
                %pdf_y = all_data.(groupName).(mouseName).(cellID).ISI_pdf_y;
                %pdf_x = all_data.(groupName).(mouseName).(cellID).ISI_pdf_x;
                %pdf_y = pdf_y(pdf_x <= 100);
                %ISI_PDFs{end+1,1} = pdf_y;

                %ISI_PDFs_X{end+1,1} = pdf_x(pdf_x <= 100);

                % if length(pdf_x) > length(x)
                %     x = pdf_x;
                % end
            end
        end
    end
end

%% plotting
figure;

% g = gramm('x',ISI_PDFs_X, 'y',ISI_PDFs, 'color',groupsVec);
% g.stat_summary('type','sem');

g = gramm('x',ISIs, 'color',groupsVec, 'column',cellTypesVec);
g.stat_density('npoints',1000, 'extra_x',0);
g.set_names('x','Inter-Spike Interval (ms)', 'Color','', 'Column','');
g.draw();

%g.axe_property('XTickLabel',xticklabels);
for ii = 1:length(unique(cellTypesVec))
    % update x tick labels for log scale
    xticklabels = g.facet_axes_handles(1,ii).XTickLabel;
    for jj = 1:length(xticklabels)
        label_jj_log = xticklabels{jj,1};
        label_jj_lin = num2str(10^str2double(label_jj_log));
        xticklabels{jj,1} = label_jj_lin;
    end
    g.facet_axes_handles(1,ii).XTickLabel = xticklabels;
end

end


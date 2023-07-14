function ISI_PDFs_compare(all_data, cell_types)
%ISI_PDFS_COMPARE Summary of this function goes here
%   Detailed explanation goes here

groupNames = fieldnames(all_data);

groupsVec = {};
ISIs = {}; % [units x times]
x = [];
ISI_PDFs_X = {};

for groupNum = 1:length(groupNames)
    groupName = groupNames{groupNum};

    mouseNames = fieldnames(all_data.(groupName));

    for mouseNum = 1:length(mouseNames)
        mouseName = mouseNames{mouseNum};

        cellIDs = fieldnames(all_data.(groupName).(mouseName));

        for cellID_num = 1:length(cellIDs)
            cellID = cellIDs{cellID_num};

            thisCellType = all_data.(groupName).(mouseName).(cellID).Cell_Type;
            isSingleUnit = all_data.(groupName).(mouseName).(cellID).IsSingleUnit;
            if any(strcmp(thisCellType, cell_types)) && isSingleUnit
                groupsVec{end+1} = groupName;
                ISI_vec = all_data.(groupName).(mouseName).(cellID).ISI_baseline_vec;
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

g = gramm('x',ISIs, 'color',groupsVec);
g.stat_density('npoints',1000, 'extra_x',0);
g.set_names('x','Inter-Spike Interval (ms)', 'Color','');
g.draw();

% update x tick labels for log scale
xticklabels = g.facet_axes_handles.XTickLabel;
for ii = 1:length(xticklabels)
    label_ii_log = xticklabels{ii,1};
    label_ii_lin = num2str(10^str2double(label_ii_log));
    xticklabels{ii,1} = label_ii_lin;
end
%g.axe_property('XTickLabel',xticklabels);
g.facet_axes_handles.XTickLabel = xticklabels;

end


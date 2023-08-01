function makeNice(axes)
% Tim C Whalen June 2019
% Make the current figure look nice for illustrator
set(axes,'box','off',...
    'TickDir','out',...
    'Fontname','Arial',... % will need to font switch to arial in AI
    'FontSize',13',...
    'LineWidth',2)
end
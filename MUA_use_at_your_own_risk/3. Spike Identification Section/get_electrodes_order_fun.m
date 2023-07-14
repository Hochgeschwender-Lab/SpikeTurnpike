function new_electrode_number = get_electrodes_order_fun(electrodes_order, n_electrodes, file_directory)

%{
    This file is for reordering electrodes based on the order provided.
%}

% The next two lines are purely for just picking out particular electrode sample
new_electrode_number = find(electrodes_order==n_electrodes);

end


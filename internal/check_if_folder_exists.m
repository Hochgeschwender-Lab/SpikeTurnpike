function check_if_folder_exists(path_directory)

%{
    This file is for checking if a particular path directory exists or not.
    If it does not exist, it creates the folder.
    Get the "path_directory" as input from the user.
%}

%% Check if folder exists, if not the "mkdir" command creates the particular "path_directory"

if not(isfolder(path_directory))
    mkdir(path_directory)
end

end
function check_if_folder_exists(folder_name)

%{
    This file is for checking if a particular named folder exists or not. If 
    it does not exist, it creates the folder.
%}

if not(isfolder(folder_name))
    mkdir(folder_name)
end

end


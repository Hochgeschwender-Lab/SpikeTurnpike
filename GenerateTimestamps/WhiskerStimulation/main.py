# Import necessary libraries
'''
os: access files stored within your computer
numpy: to perform array operations
matplotlib: to plot graphs
glob: to select only ".txt" files to analyze
scipy: to save the extracted data as ".mat: files
tkinter: to pop open the GUI to access directory
'''

import os
import matplotlib.pyplot as plt
import glob

from scipy.io import savemat
from tkinter import *
from tkinter import filedialog

# Import the necessary dependent classes defined
import data_import as di

# Clear terminal output if there are preexisting printed output
os.system('clear')

def get_file_directory_fun():
    directory_object = Tk()
    directory_object.withdraw()
    entire_file_directory = filedialog.askdirectory()
    return entire_file_directory


def find_all_folders_fun(user_selected_directory_fun):
    all_folders_names_fun = next(os.walk(user_selected_directory_fun))[1]
    
    return all_folders_names_fun


def find_all_txt_files_fun(user_selected_directory_fun):
    # Experiment file name
    find_for_txt_files = user_selected_directory_fun + "/*.txt"
    txt_file_full_paths_fun = glob.glob(find_for_txt_files)

    total_txt_files_fun = len(txt_file_full_paths_fun)
    print("There are %d number of .txt files to extract data from this folder." %total_txt_files_fun)
    
    file_names_fun = []
    
    for files in os.listdir(user_selected_directory_fun):
        # Check if the given path is a file
        if os.path.isfile(os.path.join(user_selected_directory_fun, files)):
            # Check if the file ends with a ".txt" extension
            if files.endswith(".ns6"):
                file_names_fun.append(files)
    
    return total_txt_files_fun, file_names_fun, txt_file_full_paths_fun


def extract_timestamp_info_fun(particular_txt_file_path_fun):
    print("\nExtracting data from file path: %s" %particular_txt_file_path_fun)
    session = di.Session(particular_txt_file_path_fun, int_subject_IDs=False)
    
    session_events_fun = session.events
    session_times_fun = session.times
    session_print_lines_fun = session.print_lines
    
    return session_events_fun, session_times_fun, session_print_lines_fun


def save_variables_as_mat_fun(file_directory_to_save_data_fun, variable_fun, variable_name_fun):
    print("Saving file to %s" %file_directory_to_save_data_fun)
    savemat(file_directory_to_save_data_fun + variable_name_fun + '.mat', {'data': variable_fun})


def check_if_directory_exists(directory_fun):
    if os.path.isdir(directory_fun):
        print("%s directory exists!" %directory_fun)
        
    else:
        os.mkdir(directory_fun)
    

def main():

    user_selected_directory = get_file_directory_fun()
    all_folder_names = find_all_folders_fun(user_selected_directory)
    total_num_folders = len(all_folder_names)
    print("There are %d folders in the directory" %total_num_folders)
    
    for n_folder in range(total_num_folders):
        folder_name = all_folder_names[n_folder]
        particular_file_directory = user_selected_directory + '/' + folder_name
        
        num_txt_files, file_names, txt_full_file_paths = find_all_txt_files_fun(particular_file_directory)

        for n_txtfile in range(num_txt_files):
            particular_txt_file_path = txt_full_file_paths[n_txtfile]
            particular_file_name = file_names[n_txtfile]
            
            file_directory_to_save_data = user_selected_directory + '/' + folder_name
            check_if_directory_exists(file_directory_to_save_data)
            
            file_directory_to_save_data = file_directory_to_save_data + "/MUA/"
            check_if_directory_exists(file_directory_to_save_data)
            
            file_directory_to_save_data = file_directory_to_save_data + particular_file_name[:-4]
            check_if_directory_exists(file_directory_to_save_data)
            
            file_directory_to_save_data = file_directory_to_save_data + '/allData/'
            check_if_directory_exists(file_directory_to_save_data)
            
            session_events, session_times, session_print_lines = extract_timestamp_info_fun(particular_txt_file_path)
            
            save_variables_as_mat_fun(file_directory_to_save_data, session_events, 'session_events')
            save_variables_as_mat_fun(file_directory_to_save_data, session_times, 'session_times')
            save_variables_as_mat_fun(file_directory_to_save_data, session_print_lines, 'session_print_lines')
            
            print("Finished saving file #%d" %(n_txtfile+1))
    
    return "Finished processs!"

# Run the main function  
print(main())


    
    
    





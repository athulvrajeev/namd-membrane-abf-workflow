#!/bin/bash

# Define the base directory where the files were created
base_dir="$(pwd)"   # run from the system directory

# Define the previous directory where the folders were created

# Move files to respective folders
for ((i=0; i<=10; i++))
do
    # Define the filename for the current configuration file
    filename="win${i}.conf"
    
    # Define the folder name for the current configuration file
    foldername="win${i}"
    
    # Move the file to the respective folder
    mv "${base_dir}/${filename}" "${base_dir}/${foldername}/"
done

echo "Files moved to folders successfully."


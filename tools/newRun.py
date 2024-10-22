#!/usr/bin/env python3

import argparse
import os
import shutil

def create_run_directory(run_name=None):
    #Define paths
    file_dirs = ["bin", "jobs", "analysis", "tools", "src/equil"]
    runs_dir = "runs"
    host = os.environ.get("NERSC_HOST")

    #Create the 'runs' directory if it doesn't exist
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)

    #Determine run directory name
    if run_name:
        run_dir = os.path.join(runs_dir, run_name)
    else:
        #Find the latest run directory number
        existing_runs = [d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d)) and d.startswith('run')]
        if existing_runs:
            #Remove the 'run' prefix and extract the numeric part
            last_run_number = max(int(run.lstrip('run')) for run in existing_runs if run.lstrip('run').isdigit())
            next_run_number = last_run_number + 1
        else:
            next_run_number = 1

        #Use the specified number of digits for the run directory
        run_dir = os.path.join(runs_dir, f"run{{:03d}}".format(next_run_number))

    try:
        #Create the new run directory
        os.makedirs(run_dir)

        #Copy all files and directories from 'bin' to the new run directory
        for folder in file_dirs:
            for entry in os.listdir(folder):
                full_entry_path = os.path.join(folder, entry)
                if os.path.isdir(full_entry_path):
                    #Copy the directory and its contents
                    shutil.copytree(full_entry_path, os.path.join(run_dir, entry))
                elif os.path.isfile(full_entry_path):
                    if (entry == 'gemx' or ".sh" in entry):
                        #Create an absolute symbolic link to these files
                        abs_bin_path = os.path.abspath(full_entry_path)
                        os.symlink(abs_bin_path, os.path.join(run_dir, entry))
                    #Let job script move code changes at run time so latest version.
                    #Move job scripts furhter below since only one is in bin.
                    elif (entry == 'codeChanges.txt' or entry == 'newRun.py'):
                        continue
                    elif (entry[:3] == 'job'):
                        if (host == 'perlmutter' and entry != 'job_local') or \
                        (host != 'perlmutter' and entry == 'job_local'):
                            shutil.copy(full_entry_path,run_dir)
                    else:
                        shutil.copy(full_entry_path, run_dir)

        print(f"Created run directory: {run_dir}")

    except Exception as e:
        print(f"Error occurred while creating the run directory: {e}")
        #Remove the run directory if it was created before the error occurred
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)

if __name__ == "__main__":
    #Set up and read arguments.
    parser = argparse.ArgumentParser(description="Create a new run directory and copy files.")
    parser.add_argument('-n', '--name', type=str, help="The name for the new run directory.")
    args = parser.parse_args()

    create_run_directory(run_name=args.name)

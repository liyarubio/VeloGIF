import subprocess
import os
from datetime import datetime
import pytz
import shutil
import glob
import sys
from config import Methods, output_path, figures_path, evals_path

def ensure_directory_exists(path):
    """Ensure the directory exists; create it if not."""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")
    else:
        print(f"Directory already exists: {path}")

def normalize_filename(name):
    """Normalize filenames by removing specific characters."""
    return name.lower().replace('velo', '').replace('_', '').replace('-', '').strip()

def is_file_exists(env_name, existing_files):
    """
    Check if a result file already exists.
    
    Returns:
        1 if the file exists,
        0 if the file does not exist.
    """
    env_name = env_name.lower()
    existing_lower = [f.lower() for f in existing_files]
    
    if env_name in existing_lower:
        print(f"File found: {env_name}")
        return 1  # File exists
    
    # Normalize names and check
    normalized_env = normalize_filename(env_name)
    normalized_files = [normalize_filename(f) for f in existing_files]

    if normalized_env in normalized_files:
        print(f"File found: {env_name}")
        return 1  # File exists

    # Check based on prefixes
    for file in existing_files:
        file_prefix = file.lower().split('.')[0]
        if env_name.startswith(file_prefix) or file_prefix.startswith(env_name):
            print(f"File found by prefix match: {file} (matched by {env_name})")
            return 1  # File exists
            
    print(f"Start running: {env_name}")
    return 0  # File does not exist

# Create required directories
print("\n=== Step 1: Ensuring required directories exist ===")
for path in [output_path, figures_path, evals_path]:
    print(f"Checking directory: {path}")
    ensure_directory_exists(path)

# Create log file path
log_file = os.path.join(output_path, "execution_log.txt")

# Remove h5ad extensions and collect existing files
print("\n=== Step 2: Detecting existing files ===")
existing_files = set([
    f.replace('.h5ad', '') for f in os.listdir(output_path) 
    if f.endswith(".h5ad") and not f.startswith(".")
])
print(f"Detected existing files: {existing_files}")

print("\n=== Step 3: Listing configured methods ===")
for method in Methods:
    print(f"  - {method}")

tz = pytz.timezone('Asia/Shanghai') # set timezone

with open(log_file, "a") as log:
    log.write(f"\n=== Execution started at: {datetime.now().astimezone(tz)} ===\n")
    log.write(f"Existing result files: {existing_files}\n\n")

    total_methods = len(Methods)

    for idx, env_name in enumerate(Methods, 1):
        env_name_lower = env_name.lower().strip()
        
        # Check if result file already exists
        file_exists = is_file_exists(env_name_lower, existing_files)
        if file_exists == 1:
            message = f"[{idx}/{total_methods}] Skipping {env_name}: Result file already exists\n"
            print(message)
            log.write(message)
            continue

        # Run the corresponding script
        print(f"\n[{idx}/{total_methods}] Starting to run script for: {env_name_lower}")

        run_python_file_cmd = f"conda run -n {env_name_lower} python run_{env_name_lower}.py"

        try:
            message = f"[{idx}/{total_methods}] Running environment: {env_name_lower}\n Start Time:{datetime.now().astimezone(tz)} \n"
            print(message)
            log.write(f"{message}\n")

            subprocess.run(run_python_file_cmd, shell=True, check=True)
            message = f"[{idx}/{total_methods}] Running environment: {env_name_lower}\n End Time:{datetime.now().astimezone(tz)} \n\n"
            print(message)

        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            error_message = f"Error occurred while processing {env_name}: {str(e)}\n"
            print(error_message)
            log.write(error_message)
            continue


# Run Get_Graph.py
print("\n=== Step 4: Running Get_Graph.py ===")
with open(log_file, "a") as log:
    log.write("\nStarting Get_Graph.py...\n")

    try:
        print("Starting Get_Graph.py...")
        subprocess.run("conda run -n scvelo python Get_Graph.py", 
                       shell=True, check=True)
        print("Get_Graph.py executed successfully!")
        log.write("Get_Graph.py executed successfully\n")

    except subprocess.CalledProcessError as e:
        error_message = f"Error occurred while running Get_Graph.py: {str(e)}\n"
        print(error_message)
        log.write(error_message)

    log.write(f"\n=== Execution ended at: {datetime.now().astimezone(tz)} ===\n")
    
print("\nExecution complete")
print(f"\n=== Execution ended at: {datetime.now().astimezone(tz)} ===\n")

# Delete unnecessary intermediate files
print("\n ===Delete unnecessary intermediate files=== ")
items_to_check = [
    "figures",
    "__pycache__",
    "ChEA",
    "ENCODE",
    "res",
    "saved",
    "latentvelo",
    "result/velovae",
    "tmp_model.cpt",
    "settings.txt",
    "training_loss.png"]
for folder in glob.glob('cellDancer_*'):
    if os.path.isdir(folder):
        shutil.rmtree(folder)
for item in items_to_check:
    if os.path.exists(item):
        if os.path.isdir(item): 
            shutil.rmtree(item)  
        elif os.path.isfile(item):  
            os.remove(item) 

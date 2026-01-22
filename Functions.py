
import os
import pandas as pd

def compile_era5_snowmelt(root_folder_path, output_csv="era5_snowmelt_data.csv"):

    all_dataframes = []
    
    if not os.path.exists(root_folder_path):
        print(f"Error: The folder '{root_folder_path}' does not exist.")
        return None

    sorted_folders = sorted(os.listdir(root_folder_path))

    for folder_name in sorted_folders:
        dir_path = os.path.join(root_folder_path, folder_name)

        if os.path.isdir(dir_path) and len(folder_name) == 2 and folder_name.isdigit():
            short_year = int(folder_name)
            if short_year < 50:
                full_year = f"20{folder_name}"
            else:
                full_year = f"19{folder_name}"
            
            file_name = f"era5_waves_{full_year}.nc"
            file_path = os.path.join(dir_path, file_name)

            if os.path.exists(file_path):
                try:
                    with xr.open_dataset(file_path) as ds:
                        df = ds.to_dataframe().reset_index()
                        df['source_year'] = full_year
                        all_dataframes.append(df)
                        
                except Exception as e:
                    print(f"Failed to read {file_name}: {e}")
            else:
                print(f"Warning: File {file_name} not found in folder {folder_name}")

    if all_dataframes:
        df_snowmelt = pd.concat(all_dataframes, ignore_index=True)
        df_snowmelt.to_csv(output_csv, index=False)
        return df_snowmelt
    else:
        print("No valid data found.")
        return pd.DataFrame()


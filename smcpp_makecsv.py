import json
import csv
import sys
import math
import os

# CONFIGURATION
GEN_TIME = float(os.getenv("GEN_TIME", 9.93)) #uses GEN_TIME from bash script, if not, default to 9.93

def process_file(filepath, writer, pop_name):
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        # Extract model parameters
        # smc++ JSON structure: model -> {N0, y, knots}
        model = data['model']
        N0 = model['N0']
        y_values = model['y']       # log(Ne/N0)
        knots = model['knots']      # Time points
        
        # Determine label based on filename
        filename = os.path.basename(filepath)
        if "main" in filename:
            label_type = "Main"
            unique_id = "main"
        else:
            label_type = "Bootstrap"
            # distinct ID for plotting groups (e.g. rep_0, rep_1)
            unique_id = filename.replace(".json", "").replace("_fixed", "")

        # Calculate real values and write rows
        # Time = knot * 2 * N0 * gen_time
        # Ne = exp(y) * N0
        for t_raw, y_raw in zip(knots, y_values):
            real_time = t_raw * 2 * N0 * GEN_TIME
            real_ne = math.exp(y_raw) * N0
            
            writer.writerow([pop_name, unique_id, label_type, real_time, real_ne])
            
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert_to_csv.py <output_filename.csv> <population_name> <json_files...>")
        sys.exit(1)

    output_csv = sys.argv[1]
    pop_name = sys.argv[2]
    input_files = sys.argv[3:]

    print(f"Processing {len(input_files)} files for population '{pop_name}'...")

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        # Header for R
        writer.writerow(['population', 'label', 'type', 'x', 'y'])
        
        for json_file in input_files:
            process_file(json_file, writer, pop_name)
    
    print(f"Done! Saved to {output_csv}")
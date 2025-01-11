import pandas as pd

# Columns for the dataframe
columns = ["timestep", "id", "type", "x", "y", "z", "vx", "vy", "vz"]


# Function to parse the trajectory file
def parse_trajectory(file_path):
    data = []
    timestep = None
    with open(file_path, "r") as file:
        for line in file:
            # print(line)
            if line.startswith("ITEM: TIMESTEP"):
                l = next(file).strip()
                # print(f"timestep: {l}")
                timestep = int(l)
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                num_atoms = int(next(file).strip())
                # print(f"num_atoms: {num_atoms}")
            elif line.startswith("ITEM: ATOMS"):
                for i in range(num_atoms):
                    atom_line = next(file).strip()
                    # print(f"atom_line: {atom_line}")
                    atom_data = list(map(float, atom_line.split()))
                    data.append([timestep] + atom_data)

    return pd.DataFrame(data, columns=columns)

if __name__ == "__main__":
    # Parse the provided trajectory file
    file_path = "traj.out"
    traj_df = parse_trajectory(file_path)
    # Sort by atomic ID
    sorted_traj_df = traj_df.sort_values(by=["id", "timestep"]).reset_index(drop=True)
    # Display sorted DataFrame
    sorted_traj_df.head()

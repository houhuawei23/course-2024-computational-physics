import matplotlib.pyplot as plt
import pandas as pd

# Load the uploaded file to inspect its content and structure.
file_path = "./velocity.dat"

# Read the file to determine its structure
with open(file_path, "r") as file:
    content_preview = file.readlines()[:10]  # Read the first 10 lines for a preview

content_preview

# Skip the header lines and read the data into a DataFrame
column_names = ["Chunk", "Coord1", "Ncount", "vx"]
data = pd.read_csv(file_path, delim_whitespace=True, skiprows=3, names=column_names)

# Plot vx against Coord1
plt.figure(figsize=(10, 6))
plt.plot(data["Coord1"], data["vx"], label="vx vs. Coord1", marker="o")
plt.title("Velocity (vx) vs. Coordinate (Coord1)")
plt.xlabel("Coord1")
plt.ylabel("vx")
plt.grid(True)
plt.legend()
plt.show()

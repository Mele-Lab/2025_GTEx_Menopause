import json
import re
import matplotlib.pyplot as plt

# File path
log_file = "/home/aelhomma/dino_train_logs_1000.txt"

# Initialize lists
epochs = []
train_losses = []

# Read and parse the log file
with open(log_file, "r") as f:
    for line in f:
        try:
            # Remove double commas (",,") before parsing
            cleaned_line = re.sub(r",\s*,", ",", line.strip())  
            entry = json.loads(cleaned_line)  # Parse JSON
            
            epochs.append(entry["epoch"])
            train_losses.append(entry["train_loss"])
        except json.JSONDecodeError as e:
            print(f"Skipping malformed line: {line.strip()} - Error: {e}")

# Plotting with improved aesthetics
plt.figure(figsize=(12, 6))
plt.plot(epochs, train_losses, marker="o", linestyle="-", color="b", markersize=8, linewidth=2.5)

# Titles and labels with adjusted font sizes
plt.title("Training Loss Over Epochs", fontsize=18, fontweight="bold")
plt.xlabel("Epoch", fontsize=16)
plt.ylabel("Training Loss", fontsize=16)

# Remove grid
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Save the figure before displaying
plt.savefig("convergence_plot.svg", format="svg", dpi=300)
plt.show()

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D  # Import Line2D to customize legend points
import yaml
import seaborn as sns

'''
This script does a umap for all big tissues or all tissues from the tile embeddings of the umap file:
umap file is the ouput of /gpfs/projects/bsc83/Projects/GTEx_v8/Allal/RNAPath/tiles_classification/all_tile_embeddings.py
Change the subtype_to_main mapping to choose which tissues are going to be plotted.
'''

input ="/home/laura/BSC/"
input2 ="/gpfs/projects/bsc83/Projects/GTEx_v8/"
output ="Laura/Paper_material/Post_review_analyses/UMAPs/"

# Path to the folder containing the UMAP embeddings
umap_file = input + '/Allal/RNAPath/tiles_classification/umap_embeddings_voting_selected_tissues_filt.csv'

# Initialize an empty list to hold all dataframes
umap_data_list = []
with open(input +"/Allal/RNAPath/resources/clusters.yaml", "r") as f:
    doc = yaml.load(f, Loader=yaml.FullLoader)
# Load the CSV file into a DataFrame
df = pd.read_csv(umap_file)
classes = doc['selected_tissues']['classes']
colors = doc['selected_tissues']['colors']
df['category'] = df['category'].map(lambda x: classes[int(x)] if pd.notnull(x) else x)
#df['category'] = df['top_one_label'].map(lambda x: classes[int(x)] if pd.notnull(x) else x)

# Get unique values from the 'category' column and print their length
print("Number of unique categories:", len(df['category'].unique()))


tissue_colors = {
    'breast': [73, 0, 146, 255],  # Magenta
    'ovary': [0, 146, 146, 255],  # Teal
    'vagina': [182, 109, 255, 255],  # Purple
    'uterus': [109, 182, 255, 255],  # Light Blue
    'ectocervix': [255, 182, 219, 255],  # Light Pink
    'endocervix': [255, 109, 182, 255],  # Pink
    'fallopian_tube': [219, 109, 0, 255]  # Orange
}

subtype_to_main = {
    'breast_adipocyte': 'breast',
    'breast_duct': 'breast',
    'breast_gynecomastoid_hyperplasia': 'breast',
    'breast_nerve': 'breast',
    'breast_stroma': 'breast',
    'breast_lobule': 'breast',
    'ovary_corpora': 'ovary',
    'ovary_cortex': 'ovary',
    'ovary_medulla': 'ovary',
    'ovary_vessels': 'ovary',
    'vagina_epithelium': 'vagina',
    'vagina_lamina_propria': 'vagina',
    'vagina_stroma': 'vagina',
    'vagina_vessels': 'vagina',
    'uterus_endometrium': 'uterus',
    'uterus_myometrium': 'uterus',
    'uterus_vessels': 'uterus',
    'ectocervix_epithelium': 'ectocervix',
    'ectocervix_stroma': 'ectocervix',
    'ectocervix_vessels': 'ectocervix',
    'ectocervix_glands': 'ectocervix',
    'endocervix_glandular_epithelium': 'endocervix',
    'endocervix_stroma': 'endocervix',
    'endocervix_vessels': 'endocervix',
    'fallopian_tube_epithelium': 'fallopian_tube',
    'fallopian_tube_adipose': 'fallopian_tube',
    'fallopian_tube_lumen': 'fallopian_tube',
    'fallopian_tube_vessels': 'fallopian_tube',
    'fallopian_tube_stroma': 'fallopian_tube',
    'fallopian_tube_smooth_muscle': 'fallopian_tube',
}

similarity_groups = {
    'vagina_epithelium': 'epithelium',
    'ectocervix_epithelium': 'epithelium',
    'fallopian_tube_epithelium': 'epithelium',
    'endocervix_glandular_epithelium': 'epithelium',
    'breast_adipocyte': 'adipose',
    'fallopian_tube_adipose': 'adipose',
    'fallopian_tube_vessels': 'vessels',
    'ectocervix_vessels': 'vessels',
    'endocervix_vessels': 'vessels',
    'ovary_vessels': 'vessels',
    'uterus_vessels': 'vessels',
    'vagina_vessels': 'vessels',
}


marker_shapes = {
    'epithelium': 'o',  # Circle
    'adipose': 's',     # Square
    'vessels': '^',     # Triangle
}

face_colors = {
    "uterus": [109, 182, 255, 255],         # #6DB6FF
    "ovary": [0, 146, 146, 255],            # #009292
    "vagina": [182, 109, 255, 255],         # #B66DFF
    "breast": [73, 0, 146, 255],  # #490092
    "endocervix": [255, 109, 182, 255],     # #FF6DB6
    "ectocervix": [255, 182, 219, 255],     # #FFB6DB
    "fallopian_tube": [219, 109, 0, 255]     # #DB6D00
}



df['similarity_group'] = df ['category'].map(similarity_groups)

df = df.dropna()


df['marker'] = df['similarity_group'].map(marker_shapes)

df['main_tissue'] = df['category'].map(subtype_to_main)

# Normalize the RGBA values by dividing by 255 dynamically
normalized_tissue_colors = {k: [x / 255 for x in v] for k, v in tissue_colors.items()}
normalized_face_colors = {k: [x / 255 for x in v] for k, v in face_colors.items()}

# Apply the normalized color mapping to df['edge_color']
df['face_color'] = df['main_tissue'].map(normalized_face_colors)




# Now you can plot the UMAP with the assigned colors
plt.figure(figsize=(15, 10))
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Define face colors and marker shapes
face_colors_dict = {  # Renamed to avoid conflict
    "uterus": [109, 182, 255, 255],         # #6DB6FF
    "ovary": [0, 146, 146, 255],            # #009292
    "vagina": [182, 109, 255, 255],         # #B66DFF
    "breast": [73, 0, 146, 255],            # #490092
    "endocervix": [255, 109, 182, 255],     # #FF6DB6
    "ectocervix": [255, 182, 219, 255],     # #FFB6DB
    "fallopian_tube": [219, 109, 0, 255]    # #DB6D00
}

marker_shapes = {
    'epithelium': 'o',  # Circle
    'adipose': 's',     # Square
    'vessels': '^',     # Triangle
}

#df = df[df['top_one_prob'] > 0.75]


# Remove all spines from the plot
for spine in plt.gca().spines.values():      
    spine.set_visible(False)


###################################################################################################################
##One plot per tissue
# Plot UMAP with inner (face) color based on similarity group, outer (edge) color based on tissue, and marker shape
# Loop over main tissues

marker_types = list(marker_shapes.keys())  # ['epithelium', 'adipose', 'vessels']
first_three_colors = [list(face_colors_dict.values())[i] for i in [0, 1, 3]]

# Normalize RGBA to 0-1 for matplotlib
marker_colors = {m: [c/255 for c in color] for m, color in zip(marker_types, first_three_colors)}
df['face_color'] = df['similarity_group'].map(lambda m: marker_colors[m] if m in marker_colors else (0.5, 0.5, 0.5))

for tissue in df['main_tissue'].unique():

    tissue_df = df[df['main_tissue'] == tissue]

    plt.figure(figsize=(15, 10))
    
    # Loop over subtypes within this tissue
    for subtype in tissue_df['category'].unique():
        subtype_data = tissue_df[tissue_df['category'] == subtype]
        
        # Loop over marker types within this subtype
        for marker in subtype_data['marker'].unique():
            subset = subtype_data[subtype_data['marker'] == marker]
            
# Handle face_color - check if it's a single float or a list/tuple
            face_colors_list = []  # Renamed to avoid conflict
            for color in subset['face_color']:
                if isinstance(color, (list, tuple)):  # If color is already a list/tuple, use it
                    face_colors_list.append(tuple(color))
                elif isinstance(color, (float, int)):  # If color is a single number, convert to RGB tuple
                    face_colors_list.append((color, color, color))  # Convert grayscale color
                else:
                    face_colors_list.append((0, 0, 0))  # Default to black if the color is neither list/tuple nor float
            
            plt.scatter(subset['UMAP1'], subset['UMAP2'],
                        facecolor=face_colors_list,  # Similarity group color for face
                        alpha=0.8,
                        edgecolor="black",  # Set edge color to black
                        linewidth=0.5,  # Remove edge border by setting linewidth to 0
                        s=300,  # Marker size
                        marker=marker,  # Dynamic marker shape
                        label=tissue)
    # Clean axes
    plt.xticks([]); plt.yticks([])
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    
    plt.xlabel("UMAP1", loc="left", size=20)
    plt.ylabel("UMAP2", loc="bottom", size=20)

    # Add arrows to the axes with correct length
    plt.annotate('', xy=(0.2, 0), xytext=(0, 0), 
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=3.5))  # X-axis arrow

    plt.annotate('', xy=(0, 0.2), xytext=(0, 0), 
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=3.5))  # Y-axis arrow

    plt.title(f"UMAP - {tissue.capitalize()}", fontsize=20)
    
    # Legends
    # Create custom legend handles for face colors (similarity groups) - dots with fill

    face_color_legend = [
        Line2D(
            [0], [0],
            marker=marker_shapes[m],                # marker shape
            color='w',                              # line color (ignored)
            markerfacecolor=tuple(marker_colors[m]), # convert list -> tuple
            markeredgecolor='b',                # optional border
            markersize=10,
            label=m
        )
        for m in marker_colors.keys()
    ]

    # # Add the legend
    # legend1 = plt.legend(
    #     handles=face_color_legend,
    #     title="Marker Type",
    #     loc='upper right',
    #     bbox_to_anchor=(1.3, 0.8),
    #     prop={'size': 14}  # adjust font size
    # )

    # plt.gca().add_artist(legend1)
    plt.legend(handles=face_color_legend, loc='upper right', fontsize=14, title="")

    # Add labels and title to match the previous plot
    plt.xlabel("UMAP1", loc="left", size=40)
    plt.ylabel("UMAP2", loc="bottom", size=40)

    # Add arrows to the axes with correct length (same as original plot)
    plt.annotate('', xy=(0.2, 0), xytext=(0, 0), 
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=2.5))  # X-axis arrow

    plt.annotate('', xy=(0, 0.2), xytext=(0, 0), 
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=2.5))  # Y-axis arrow
    # Final adjustments
    plt.tight_layout()
    # Save each plot
    out_path = os.path.join(input,output, f"umap_{tissue}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight", format="png")
    plt.close()

print("One plot per tissue has been saved!")












###################################################################################################################+
#Plot in main
# Plot UMAP with inner (face) color based on similarity group, outer (edge) color based on tissue, and marker shape
for subtype in df['category'].unique():
    tissue_data = df[df['category'] == subtype]
    
    # Plot with edge and face colors, and dynamic marker shapes
    for marker in tissue_data['marker'].unique():  # Iterate over unique marker shapes
        subset = tissue_data[tissue_data['marker'] == marker]
        
        # Handle face_color - check if it's a single float or a list/tuple
        face_colors_list = []  # Renamed to avoid conflict
        for color in subset['face_color']:
            if isinstance(color, (list, tuple)):  # If color is already a list/tuple, use it
                face_colors_list.append(tuple(color))
            elif isinstance(color, (float, int)):  # If color is a single number, convert to RGB tuple
                face_colors_list.append((color, color, color))  # Convert grayscale color
            else:
                face_colors_list.append((0, 0, 0))  # Default to black if the color is neither list/tuple nor float
        
        plt.scatter(subset['UMAP1'], subset['UMAP2'],
                    facecolor=face_colors_list,  # Similarity group color for face
                    alpha=0.8,
                    edgecolor="black",  # Set edge color to black
                    linewidth=0.5,  # Remove edge border by setting linewidth to 0
                    s=70,  # Marker size
                    marker=marker,  # Dynamic marker shape
                    label=subtype)

# Remove ticks from the axes
plt.xticks([])  # Remove x-axis tick labels
plt.yticks([])  # Remove y-axis tick labels
# Remove specific spines (e.g., top and right)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)

# Add arrows to the axes with correct length
plt.annotate('', xy=(0.2, 0), xytext=(0, 0), 
             xycoords='axes fraction', textcoords='axes fraction',
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=3.5))  # X-axis arrow

plt.annotate('', xy=(0, 0.2), xytext=(0, 0), 
             xycoords='axes fraction', textcoords='axes fraction',
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=3.5))  # Y-axis arrow

# Create custom legend handles for face colors (similarity groups) - dots with fill
face_color_legend = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=[c/255 for c in color], 
           markersize=10, label=name) for name, color in face_colors_dict.items()
]

# Create custom legend handles for edge colors (tissues) - dots with no fill, only border
edge_color_legend = [
    plt.Line2D([0], [0], marker=shape, color='w', markerfacecolor="none", markeredgecolor="black", 
               markersize=8, label=label)  # Smaller marker size
    for label, shape in marker_shapes.items()
]


# Create the first legend (Face color: Similarity Groups)
legend1 = plt.legend(handles=face_color_legend, title="", loc='upper right', 
                     bbox_to_anchor=(1.3, 0.8), prop={'size': 20})  # Font size 20, same as original plot

# Create the second legend (Edge color: Tissues) at a different vertical position
legend2 = plt.legend(handles=edge_color_legend, title="", loc='upper right', 
                     bbox_to_anchor=(1.3, 1), prop={'size': 20})  # Font size 20, same as original plot

# Add the first legend to the plot
plt.gca().add_artist(legend1)

# Add labels and title to match the previous plot
plt.xlabel("UMAP1", loc="left", size=20)
plt.ylabel("UMAP2", loc="bottom", size=20)

# Add arrows to the axes with correct length (same as original plot)
plt.annotate('', xy=(0.2, 0), xytext=(0, 0), 
             xycoords='axes fraction', textcoords='axes fraction',
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=2.5))  # X-axis arrow

plt.annotate('', xy=(0, 0.2), xytext=(0, 0), 
             xycoords='axes fraction', textcoords='axes fraction',
             arrowprops=dict(facecolor='black', arrowstyle='->', lw=2.5))  # Y-axis arrow

# Final adjustments
plt.tight_layout()
plt.show() 
plt.savefig(input+"/Laura/Paper_material/Post_review_analyses/4umap_main_fig_probs.png", format="png", dpi=300)
plt.close()  # Close the plot to free up memory



plt.figure(figsize=(10, 8))

# Scatter plot of UMAP1 vs. UMAP2, color by 'top_one_prob'
scatter = plt.scatter(df['UMAP1'], df['UMAP2'], c=df['top_one_prob'], cmap='viridis', s=10, alpha=0.7)

# Add color bar
plt.colorbar(scatter, label='Top One Probability')

# Set axis labels and title
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP Colored by Top One Probability')

# Show the plot
plt.savefig("/gpfs/projects/bsc83/Projects/GTEx_v8/Allal/RNAPath/tiles_classification/umap_top1_prob_4.png")



plt.figure(figsize=(10, 8))

# Scatter plot of UMAP1 vs. UMAP2, color by 'top_one_prob'
scatter = plt.scatter(df['UMAP1'], df['UMAP2'], c=df['top_two_prob'], cmap='viridis', s=10, alpha=0.7)

# Add color bar
plt.colorbar(scatter, label='Top Two Probability')

# Set axis labels and title
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP Colored by Top Two Probability')

# Show the plot
plt.savefig("/gpfs/projects/bsc83/Projects/GTEx_v8/Allal/RNAPath/tiles_classification/umap_top2_prob_4.png")




import matplotlib.pyplot as plt

# Filter points with top_one_prob > 0.8
filtered_df = df[df['top_one_prob'] > 0.8]

plt.figure(figsize=(10, 8))

# Scatter plot of filtered UMAP points
scatter = plt.scatter(
    filtered_df['UMAP1'], 
    filtered_df['UMAP2'], 
    c=filtered_df['top_one_prob'], 
    cmap='viridis', 
    s=10, 
    alpha=0.7
)

# Add color bar
plt.colorbar(scatter, label='Top One Probability')

# Set axis labels and title
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP Colored by Top One Probability (top_one_prob > 0.8)')

# Save plot
plt.savefig("/gpfs/projects/bsc83/Projects/GTEx_v8/Allal/RNAPath/tiles_classification/umap_top1_prob_filtered_4.png")

print("Filtered UMAP plot saved successfully!")




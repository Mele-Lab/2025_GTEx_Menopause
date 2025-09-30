import os
import glob
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

# 1. Folder containing all your CSVs
folder = '/X/files/Ovary_ALL750'

# 2. Prepare the “true 8‐neighborhood” sign‐pairs
expected_offsets = {
    (dx, dy)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    if not (dx == 0 and dy == 0)
}

# 3. Prepare containers
file_counts = {}  # map filename -> number of interior tiles
counts_list = []  # list of counts in the same order

# 4. Iterate CSVs
for path in sorted(glob.glob(os.path.join(folder, '*.csv'))):
    # Load & reset index
    df = pd.read_csv(path).reset_index(drop=True)
    coords = df[['coord_x', 'coord_y']].values

    # Build NN model and get indices
    nbrs = NearestNeighbors(n_neighbors=9, algorithm='auto').fit(coords)
    _, indices = nbrs.kneighbors(coords)

    # 4a. Select candidates: label==4 & neighbor‐labels ⊆ {2,4}
    candidates = []
    for i in df.index[df['label'] == 4]:
        neigh_idxs = indices[i][1:9]
        neigh_labels = df['label'].iloc[neigh_idxs].tolist()
        if set(neigh_labels).issubset({2, 4}):
            candidates.append(i)

    # 4b. Filter out edge tiles by inferring grid spacing “a”
    interior = []
    for i in candidates:
        x0, y0 = coords[i]
        neigh_coords = coords[indices[i][1:9]]
        deltas = neigh_coords - np.array([x0, y0])

        # infer spacing from axis‐aligned neighbors
        axis_aligned = []
        for dx, dy in deltas:
            if abs(dx) < 1e-6 and abs(dy) > 1e-6:
                axis_aligned.append(abs(dy))
            elif abs(dy) < 1e-6 and abs(dx) > 1e-6:
                axis_aligned.append(abs(dx))
        if not axis_aligned:
            continue

        a = float(np.median(axis_aligned))
        # normalize & round
        norm = {
            (int(round(dx / a)), int(round(dy / a)))
            for dx, dy in deltas
        }

        if norm == expected_offsets:
            interior.append(i)

    # 4c. Record results
    fname = os.path.basename(path)
    count = len(interior)
    file_counts[fname] = count
    counts_list.append(count)

# 5. Output
print("Counts of true‐interior tiles per file:")
for fname, cnt in file_counts.items():
    print(f"  {fname}: {cnt}")
print("\nAs list of counts:", counts_list)

import os, glob
import pandas as pd

# 1. Load your metadata (only need Subject.ID and Age.Bracket)
meta = pd.read_csv('/home/dtabares/Escritorio/files/histological_data.csv', usecols=['Subject.ID', 'Age.Bracket'])

# Build a lookup: Subject.ID → Age.Bracket
# If there are multiple rows per subject, we'll just take the first
subject_to_bracket = meta.drop_duplicates('Subject.ID').set_index('Subject.ID')['Age.Bracket'].to_dict()

# Define how to encode each bracket as an integer
age_map = {
    '20-29': 0,
    '30-39': 1,
    '40-49': 2,
    '50-59': 3,
    '60-69': 4,
    '70-79': 5
}

# 2. Gather your slide CSVs
folder    = '/home/dtabares/Escritorio/files/Ovary_ALL750'
csv_paths = sorted(glob.glob(os.path.join(folder, '*.csv')))

# 3. Build the list of age‑group labels
age_labels = []
for path in csv_paths:
    fname = os.path.basename(path).rsplit('.', 1)[0]       # e.g. "GTEX-1AMFI-1425"
    subj  = '-'.join(fname.split('-', 2)[:2])              # e.g. "GTEX-1AMFI"
    
    bracket = subject_to_bracket.get(subj)
    if bracket is None:
        # no match in metadata
        age_labels.append(None)
    else:
        # map "60-69" → 4, etc.
        age_labels.append(age_map.get(bracket, None))

# 4. Inspect your result
print(age_labels)

for i, (path, label) in enumerate(zip(csv_paths, age_labels)):
    fname = os.path.basename(path)
    print(f"{i}. {fname} → age‐group code: {label}")



import matplotlib.pyplot as plt
import pandas as pd
import os
# 1. Build DataFrame
df_plot = pd.DataFrame({
    'slide_idx': range(len(counts_list)),
    'count':      counts_list,
    'age_group':  age_labels
}).dropna(subset=['age_group'])
df_plot['age_group'] = df_plot['age_group'].astype(int)

# 2. Map codes back to bracket strings
code_to_bracket = {
    0: '20-29',
    1: '30-39',
    2: '40-49',
    3: '50-59',
    4: '60-69',
    5: '70-79'
}

# 3. Plot
cmap = plt.cm.get_cmap('tab10', df_plot['age_group'].nunique())
fig, ax = plt.subplots(figsize=(8,5))

for code in sorted(df_plot['age_group'].unique()):
    subset = df_plot[df_plot['age_group'] == code]
    ax.scatter(
        subset['slide_idx'],
        subset['count'],
        s=30,                       
        color=cmap(code),
        label=code_to_bracket[code]
    )

ax.set_xlabel('Slide index')
ax.set_ylabel('Filtered number of tiles of follicles')
ax.legend(title='Age group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()




import pandas as pd
import matplotlib.pyplot as plt

# 1. Build DataFrame
df = pd.DataFrame({
    'count':      counts_list,
    'age_group':  age_labels
}).dropna(subset=['age_group'])
df['age_group'] = df['age_group'].astype(int)

# 2. Map integer codes back to bracket labels
code_to_bracket = {
    0: '20-29',
    1: '30-39',
    2: '40-49',
    3: '50-59',
    4: '60-69',
    5: '70-79'
}
df['bracket'] = df['age_group'].map(code_to_bracket)

# 3. Aggregate total counts and sample sizes per bracket
agg_counts = df.groupby('bracket')['count'].sum()
sample_sizes = df.groupby('bracket').size()

# Ensure we follow the defined bracket order
ordered_brackets = list(code_to_bracket.values())
agg_counts = agg_counts.reindex(ordered_brackets, fill_value=0)
sample_sizes = sample_sizes.reindex(ordered_brackets, fill_value=0)

# 4. Build x‐labels including sample size
x_labels = [f"{br} (n={sample_sizes[br]})" for br in ordered_brackets]

# 5. Plot
fig, ax = plt.subplots(figsize=(8,5))
ax.bar(x_labels, agg_counts.values)
ax.set_xlabel('Age group')
ax.set_ylabel('Total number of tiles of follicles')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.show()




import pandas as pd
import matplotlib.pyplot as plt

# 1. Build your DataFrame (reuse counts_list & age_labels)
df = pd.DataFrame({
    'count':     counts_list,
    'age_group': age_labels
}).dropna(subset=['age_group'])
df['age_group'] = df['age_group'].astype(int)

# 2. Map codes → bracket strings
code_to_bracket = {
    0: '20-29',
    1: '30-39',
    2: '40-49',
    3: '50-59',
    4: '60-69',
    5: '70-79'
}
df['bracket'] = df['age_group'].map(code_to_bracket)

# 3. Compute mean interior-tile count per bracket
#    (this is total count / number of samples in each bracket)
means = df.groupby('bracket')['count'].mean()

# 4. Reindex to ensure proper age order
ordered = ['20-29','30-39','40-49','50-59','60-69','70-79']
means = means.reindex(ordered)

# 5. Plot bar chart
fig, ax = plt.subplots(figsize=(8,5))
bars = ax.bar(means.index, means.values)

# 6. Annotate each bar with its mean value
for bar in bars:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width()/2,
        h + 2,                    # small offset above bar
        f"{h:.1f}",
        ha='center',
        va='bottom'
    )

ax.set_xlabel('Age group')
ax.set_ylabel('Mean number of tiles of follicles')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.show()



import pandas as pd
import matplotlib.pyplot as plt

# --- assume counts_list and age_labels already exist ---

# 1. Build DataFrame
df = pd.DataFrame({
    'count':      counts_list,
    'age_group':  age_labels
}).dropna(subset=['age_group'])
df['age_group'] = df['age_group'].astype(int)

# 2. Map integer codes → bracket strings
code_to_bracket = {
    0: '20-29',
    1: '30-39',
    2: '40-49',
    3: '50-59',
    4: '60-69',
    5: '70-79'
}
df['bracket'] = df['age_group'].map(code_to_bracket)

# 3. Compute median interior‑tile count per bracket
ordered = ['20-29','30-39','40-49','50-59','60-69','70-79']
medians = df.groupby('bracket')['count'].median().reindex(ordered)

# 4. Plot bar chart of medians
fig, ax = plt.subplots(figsize=(8,5))
bars = ax.bar(medians.index, medians.values)

# 5. Annotate each bar with its median value
for bar in bars:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width()/2,
        h + max(medians.values)*0.01,   # small offset above bar
        f"{h:.1f}",
        ha='center',
        va='bottom'
    )

ax.set_xlabel('Age groups')
ax.set_ylabel('Median number of tiles of follicles')
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
plt.show()



import pandas as pd
import matplotlib.pyplot as plt

# 1. Build DataFrame from your existing lists
df = pd.DataFrame({
    'count':      counts_list,
    'age_group':  age_labels
}).dropna(subset=['age_group'])
df['age_group'] = df['age_group'].astype(int)

# 2. Map integer codes → original bracket strings
code_to_bracket = {
    0: '20-29',
    1: '30-39',
    2: '40-49',
    3: '50-59',
    4: '60-69',
    5: '70-79'
}
df['bracket'] = df['age_group'].map(code_to_bracket)

# 3. Define your 3 larger age‐ranges
def super_group(br):
    if br in ('20-29', '30-39'):
        return '20-39'
    if br in ('40-49', '50-59'):
        return '40-59'
    return '60-79'

df['group'] = df['bracket'].apply(super_group)

# 4. Compute mean interior‐tile count per super‐group
order = ['20-39', '40-59', '60-79']
means = df.groupby('group')['count'].mean().reindex(order)

# 5. Plot
fig, ax = plt.subplots(figsize=(6,4))
bars = ax.bar(means.index, means.values)

# 6. Annotate each bar with its mean value
for bar in bars:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width()/2,
        h + h*0.02,          # 2% above bar
        f"{h:.1f}",
        ha='center',
        va='bottom'
    )

ax.set_xlabel('Age groups')
ax.set_ylabel('Mean number of tiles of follicles')
plt.tight_layout()
plt.show()





import os
import glob
import pandas as pd

def main():
    # 1. Path to the GTEx phenotype file (gzipped)
    phen_file = '/home/dtabares/Escritorio/files/GTEx_Subject_Phenotypes.GRU.txt.gz'

    # 2. Load GTEx phenotypes and build lookup: SUBJID -> AGE
    phen_df = pd.read_csv(
        phen_file,
        sep='\t',
        compression='gzip',
        low_memory=False
    )
    phen_map = (phen_df
                .drop_duplicates('SUBJID')
                .set_index('SUBJID')['AGE']
                .to_dict())

    # 3. Folder containing your slide files
    folder = '/home/dtabares/Escritorio/files/Ovary_ALL750'

    # 4. Gather all CSV files and extract donor IDs
    csv_paths = sorted(glob.glob(os.path.join(folder, '*.csv')))
    donor_ids = []
    for path in csv_paths:
        fname = os.path.basename(path).rsplit('.', 1)[0]  # e.g. "GTEX-1AMFI-1425"
        donor = '-'.join(fname.split('-', 2)[:2])         # e.g. "GTEX-1AMFI"
        donor_ids.append(donor)

    # 5. Map each donor ID to its age (None if not found)
    ages = [phen_map.get(donor) for donor in donor_ids]

    # 6. Create DataFrame of results
    results_df = pd.DataFrame({
        'donor_id': donor_ids,
        'age':      ages
    })

    # 7. Output the results
    print(results_df.to_string(index=False))

    # Optionally, save to CSV:
    # results_df.to_csv('donor_ages.csv', index=False)
    print(ages)

if __name__ == '__main__':
    main()

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

# 1. Folder containing all your CSVs
folder = '/home/dtabares/Escritorio/files/Ovary_ALL750'

# 2. Prepare the “true 8‑neighborhood” sign‑pairs
expected_offsets = {
    (dx, dy)
    for dx in (-1, 0, 1)
    for dy in (-1, 0, 1)
    if not (dx == 0 and dy == 0)
}

# 3. Containers
counts_list = []
ages = []

# 4. Load GTEx phenotypes once
phen_file = '/home/dtabares/Escritorio/files/GTEx_Subject_Phenotypes.GRU.txt.gz'
phen = pd.read_csv(phen_file, sep='\t', compression='gzip', low_memory=False)
phen_map = phen.drop_duplicates('SUBJID').set_index('SUBJID')['AGE'].to_dict()

# 5. Iterate each slide CSV
for path in sorted(glob.glob(os.path.join(folder, '*.csv'))):
    # --- record age ---
    fname = os.path.basename(path).rsplit('.',1)[0]   # "GTEX-XXXX-YYYY"
    donor = '-'.join(fname.split('-',2)[:2])          # "GTEX-XXXX"
    age = phen_map.get(donor, np.nan)
    ages.append(age)

    # --- compute interior‑tile count ---
    df = pd.read_csv(path).reset_index(drop=True)
    coords = df[['coord_x','coord_y']].values
    nbrs = NearestNeighbors(n_neighbors=9).fit(coords)
    _, idxs = nbrs.kneighbors(coords)

    # pick label==4 candidates whose 8-NN labels ⊆ {2,4}
    cands = [i for i in df.index[df['label']==4]
             if set(df['label'].iloc[idxs[i][1:9]].tolist()).issubset({2,4})]

    interior = []
    for i in cands:
        x0,y0 = coords[i]
        deltas = coords[idxs[i][1:9]] - np.array([x0,y0])
        # infer spacing
        axis = [abs(dy) if abs(dx)<1e-6 else abs(dx)
                for dx,dy in deltas if (abs(dx)<1e-6) ^ (abs(dy)<1e-6)]
        if not axis:
            continue
        a = np.median(axis)
        norm = {(int(round(dx/a)), int(round(dy/a))) for dx,dy in deltas}
        if norm == expected_offsets:
            interior.append(i)

    counts_list.append(len(interior))

# 6. Build DataFrame & drop missing ages
df_plot = pd.DataFrame({'age': ages, 'count': counts_list}).dropna()
df_plot['age'] = df_plot['age'].astype(int)

# 7. Plot count vs age
plt.figure(figsize=(8,5))
plt.scatter(df_plot['age'], df_plot['count'], s=10, alpha=0.7)
plt.xlabel('Donor Age (years)')
plt.ylabel('Total number of tiles of follicles')
plt.tight_layout()
plt.show()



import pandas as pd
import matplotlib.pyplot as plt

# -- assume you have these two lists from your earlier loop --
# ages        = [...]  # one age (int) or np.nan per sample
# counts_list = [...]  # one interior‑tile count per sample

# 1. Build DataFrame and drop missing ages
df = pd.DataFrame({
    'age':  ages,
    'count': counts_list
}).dropna(subset=['age'])
df['age'] = df['age'].astype(int)

# 2. Compute mean count per age
means = df.groupby('age')['count'].mean().sort_index()

# 3. Plot
fig, ax = plt.subplots(figsize=(8,5))
ax.bar(means.index.astype(str), means.values)
ax.set_xlabel('Donor age (years)')
ax.set_ylabel('Mean number of tiles of follicles')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


import pandas as pd
import matplotlib.pyplot as plt

# === Immediately after your main loop, where ages & counts_list are populated ===

# 1. Assemble DataFrame
df = pd.DataFrame({
    'age':  ages,
    'count': counts_list
}).dropna(subset=['age'])
df['age'] = df['age'].astype(int)

# 2. Median counts per exact age
medians = df.groupby('age')['count'].median().sort_index()

# 3. Plot bars only
fig, ax = plt.subplots(figsize=(8,5))
ax.bar(medians.index.astype(str), medians.values)

ax.set_xlabel('Donor age (years)')
ax.set_ylabel('Median number of tiles of follicles')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- 1) Count cortex tiles (label == 2) per slide ---
cortex_counts = []
for path in csv_paths:
    df_slide = pd.read_csv(path, usecols=['label'])
    cortex_counts.append(int((df_slide['label'] == 2).sum()))

# --- 2) Build dataframe with exact ages ---
df = pd.DataFrame({'age': ages, 'cortex': cortex_counts}).dropna(subset=['age'])
df['age'] = df['age'].astype(int)

# --- 3) Bin ages into 10-year groups ---
bins   = [20, 30, 40, 50, 60, 70, 80]              # [20–29], [30–39], ...
labels = ['20-29','30-39','40-49','50-59','60-69','70-79']
df['decade'] = pd.cut(df['age'], bins=bins, right=False, labels=labels)
df = df.dropna(subset=['decade'])

# Collect data in plotting order
order = labels
data  = [df.loc[df['decade'] == br, 'cortex'].values for br in order]

# --- 4) Boxplot + jittered points (same styling as before) ---
fig, ax = plt.subplots(figsize=(6,6))
bp = ax.boxplot(
    data, labels=order, patch_artist=True,
    showfliers=False, whis=1.5
)

for box in bp['boxes']:
    box.set(facecolor='#C084FC80', edgecolor='black', linewidth=1.5)
for median in bp['medians']:
    median.set(color='black', linewidth=2)
for whisker in bp['whiskers']:
    whisker.set(color='black', linewidth=1.2)
for cap in bp['caps']:
    cap.set(color='black', linewidth=1.2)

# jitter every slide
for i, ys in enumerate(data, start=1):
    xs = np.random.normal(loc=i, scale=0.06, size=len(ys))
    ax.scatter(xs, ys, s=18, alpha=0.7, color='#8B5CF6', edgecolors='none')

# show n per bin
ns = [len(v) for v in data]
ax.set_xticklabels([f"{br}\n(n={n})" for br, n in zip(order, ns)])

ax.set_xlabel('Age (years)')
ax.set_ylabel('Cortex tiles')
ax.set_ylim(bottom=0)  # Optional: if the spread is huge, try log scale:
# ax.set_yscale('log'); ax.set_ylabel('Cortex tiles (log scale)')
plt.tight_layout()
plt.show()




import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

folder = '/home/dtabares/Escritorio/files/Ovary_ALL750'
csv_paths = sorted(glob.glob(os.path.join(folder, '*.csv')))

cortex_counts = []
for path in csv_paths:
    df_tmp = pd.read_csv(path, usecols=['label'])      # lightweight read
    cortex_counts.append((df_tmp['label'] == 2).sum())

# sanity check: interior & cortex lists must align
assert len(cortex_counts) == len(counts_list) == len(age_labels)

# -------------------------------------------------------------------
# 2. Build the normalised metric
# -------------------------------------------------------------------
normalised_counts = [
    (fol/count) if count > 0 else np.nan      # avoid divide-by-zero
    for fol, count in zip(counts_list, cortex_counts)
]

# -------------------------------------------------------------------
# 3. Tidy DataFrame for plotting
# -------------------------------------------------------------------
df_norm = (
    pd.DataFrame({
        'norm_count':  normalised_counts,
        'age_group':   age_labels
    })
    .dropna(subset=['age_group', 'norm_count'])
    .astype({'age_group': int})
)

code_to_bracket = {
    0: '20-29', 1: '30-39', 2: '40-49',
    3: '50-59', 4: '60-69', 5: '70-79'
}
ordered = list(code_to_bracket.values())
df_norm['bracket'] = (
    df_norm['age_group']
      .map(code_to_bracket)
      .astype(pd.CategoricalDtype(ordered, ordered=True))
)

# -------------------------------------------------------------------
# 4. Box-and-scatter plot (normalised values)
# -------------------------------------------------------------------
sns.set_style('whitegrid')
fig, ax = plt.subplots(figsize=(6,4))

sns.boxplot(
    data=df_norm,
    x='bracket', y='norm_count',
    order=ordered, width=.6, fliersize=0,
    boxprops=dict(facecolor='#c9a3ff66'),
    medianprops=dict(color='k', linewidth=2),
    ax=ax
)

sns.stripplot(
    data=df_norm,
    x='bracket', y='norm_count',
    order=ordered, size=4, jitter=.25,
    alpha=.6, color='#7b42ff', ax=ax
)

ax.set_xlabel('Age (years)')
ax.set_ylabel('Interior-follicle tiles per cortex tile')
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv("species_similarity_matrix.csv", index_col=0)

# Force both index and columns to string
df.index = df.index.astype(str)
df.columns = df.columns.astype(str)

order = df.mean(axis=1).sort_values(ascending=False).index # dataframe is ordered based on the columns values
df_ordered = df.loc[order, order]

plt.figure(figsize=(10,8))
sns.heatmap(
    df_ordered.astype(float),
    cmap="viridis",
    square=True,
    cbar_kws={"label": "Similarity (%)"},
    annot=True,
    fmt=".1f",
    annot_kws={"size": 7, "color": "black"}
)
plt.title("Species Similarity Heatmap (annotated)")
plt.tight_layout()
plt.savefig("species_similarity_heatmap_annotated.png", dpi=300)

# Spatial Transcriptomic Analysis of Glioblastoma

## Overview

This project explores the spatial transcriptomic landscape of **human glioblastoma (GBM)** using a **10x Genomics Visium Spatial Gene Expression dataset**. The analysis combines gene-expression profiles, spatial coordinates, and histological tissue information to investigate transcriptomic heterogeneity across the tumor tissue.

The workflow includes:

- quality control and filtering,
- normalization and log transformation,
- highly variable gene selection,
- PCA and UMAP dimensionality reduction,
- Leiden clustering,
- spatial visualization of transcriptomic clusters,
- cluster-level differential gene expression,
- spatially variable gene detection using SpatialDE,
- and literature-guided visualization of previously reported GBM-related genes.

The project is implemented mainly with **Scanpy** and **SpatialDE**.

---

## Dataset

The dataset is loaded directly using Scanpy:

```python
adata = sc.datasets.visium_sge(
    sample_id="Parent_Visium_Human_Glioblastoma"
)
```

The original dataset contains:

| Property | Value |
|---|---:|
| Spatial spots | 3,468 |
| Genes | 36,601 |
| Technology | 10x Genomics Visium |
| Tissue | Human Glioblastoma |
| Data type | Spatial Gene Expression |

Each Visium spot contains a gene-expression profile together with its spatial location within the tissue.

The initial data matrix is therefore:

```text
3,468 spatial spots × 36,601 genes
```

In addition to the expression matrix, the AnnData object contains:

- spot metadata in `adata.obs`,
- gene metadata in `adata.var`,
- spatial coordinates in `adata.obsm["spatial"]`,
- and tissue-image information in `adata.uns["spatial"]`.

---

## Project Objectives

The main objectives are to:

1. assess and filter low-quality spatial transcriptomic observations,
2. reduce technical variability through normalization,
3. identify informative genes,
4. characterize transcriptomic heterogeneity using unsupervised analysis,
5. detect transcriptionally distinct spatial regions,
6. map transcriptomic clusters back onto the tissue image,
7. identify cluster-associated marker genes,
8. investigate biologically relevant genes such as `EGFR`,
9. detect spatially variable genes using SpatialDE,
10. examine the spatial distribution of previously reported GBM subtype-associated genes.

---

## Analysis Workflow

```text
Human Glioblastoma Visium Dataset
            │
            ▼
      Quality Control
            │
            ▼
         Filtering
            │
            ▼
Normalization + Log Transformation
            │
            ▼
 Highly Variable Gene Selection
            │
            ▼
           PCA
            │
            ▼
 Nearest-Neighbor Graph
            │
            ▼
          UMAP
            │
            ▼
   Leiden Clustering
            │
            ▼
 Spatial Cluster Mapping
            │
            ▼
 Differential Expression
            │
            ▼
      SpatialDE
            │
            ▼
Spatially Variable Genes
            │
            ▼
Literature-Guided GBM Gene Analysis
```

---

## 1. Data Loading and Gene Name Handling

The dataset is loaded as an `AnnData` object using Scanpy.

Because duplicated gene names may exist in the original dataset, gene names are made unique:

```python
adata.var_names_make_unique()
```

This prevents ambiguity during gene indexing and downstream analysis.

---

## 2. Quality Control

Mitochondrial genes are identified using the `MT-` prefix:

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
```

Quality-control metrics are then calculated:

```python
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)
```

The main QC variables include:

- `total_counts` — total transcript count per spatial spot,
- `n_genes_by_counts` — number of detected genes per spot,
- `pct_counts_mt` — percentage of mitochondrial transcripts.

The distributions of these variables are visualized before applying filtering thresholds.

---

## 3. Spot and Gene Filtering

Low-quality or extreme observations are removed using the following criteria:

```python
sc.pp.filter_cells(adata, min_counts=2000)
sc.pp.filter_cells(adata, max_counts=20000)

adata = adata[
    adata.obs["pct_counts_mt"] < 20
]

sc.pp.filter_genes(
    adata,
    min_cells=10
)
```

After filtering, the dataset contains approximately:

```text
2,870 spatial spots × 18,618 genes
```

This removes low-information spots, extreme high-count observations, spots with elevated mitochondrial expression, and genes detected in very few spatial locations.

---

## 4. Normalization and Log Transformation

Gene-expression values are normalized to reduce variability caused by differences in sequencing depth:

```python
sc.pp.normalize_total(
    adata,
    inplace=True
)
```

The normalized expression values are then log-transformed:

```python
sc.pp.log1p(adata)
```

This applies:

```text
log(1 + x)
```

to the expression values and reduces the influence of extremely highly expressed genes.

---

## 5. Highly Variable Gene Selection

Highly variable genes are selected using the Seurat-based method:

```python
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=200
)
```

These genes contain relatively strong variation across spatial spots and provide an informative representation of the transcriptomic structure.

For the main clustering workflow, the top **200 highly variable genes** are used.

---

## 6. Dimensionality Reduction

### Principal Component Analysis

PCA is used to reduce the high-dimensional gene-expression space:

```python
sc.pp.pca(adata)
```

The resulting principal components summarize the major sources of transcriptomic variation.

### Neighborhood Graph

A nearest-neighbor graph is constructed:

```python
sc.pp.neighbors(adata)
```

The graph represents similarity between spots in transcriptomic space rather than direct physical distance within the tissue.

### UMAP

The transcriptomic structure is projected into two dimensions:

```python
sc.tl.umap(adata)
```

Each point represents one Visium spot, with nearby points generally having more similar expression profiles.

---

## 7. Leiden Clustering

Unsupervised graph-based clustering is performed using the Leiden algorithm:

```python
sc.tl.leiden(
    adata,
    key_added="clusters"
)
```

The analysis identifies **11 transcriptomically distinct clusters**:

```text
0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
```

These clusters should initially be interpreted as groups of spatial spots with similar transcriptomic profiles.

They are not automatically equivalent to:

- GBM molecular subtypes,
- individual cell types,
- or histopathological regions.

Further biological annotation is required.

---

## 8. UMAP Visualization

The UMAP structure is examined using:

```python
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "clusters"
    ],
    wspace=0.4
)
```

This allows comparison between cluster structure and technical quality-control variables such as sequencing depth and gene detection.

---

## 9. Spatial Visualization

QC variables are mapped directly onto the tissue:

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color=[
        "total_counts",
        "n_genes_by_counts"
    ]
)
```

The Leiden clusters are then projected onto the histological tissue image:

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    size=1.5,
    alpha=0.8
)
```

This step connects transcriptomic heterogeneity with physical tissue organization.

Several clusters show spatially coherent distributions, suggesting that transcriptional heterogeneity is regionally structured rather than randomly distributed.

---

## 10. Differential Gene Expression

Cluster-associated marker genes are identified using:

```python
sc.tl.rank_genes_groups(
    adata,
    "clusters",
    method="t-test"
)
```

For example, the top genes associated with Cluster 1 are visualized using:

```python
sc.pl.rank_genes_groups_heatmap(
    adata,
    groups="1",
    n_genes=10,
    groupby="clusters"
)
```

Genes observed among the cluster-associated markers include examples such as:

```text
PTN
S100B
FABP5
MT3
DBI
APLP2
MT2A
METTL7B
TIMP1
```

These genes may help characterize the molecular identity of specific transcriptomic regions.

---

## 11. EGFR Spatial Expression

The spatial distribution of `EGFR`, an important gene in glioblastoma biology, is visualized alongside the discovered clusters:

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color=[
        "clusters",
        "EGFR"
    ]
)
```

The resulting map demonstrates heterogeneous EGFR expression across the tissue.

This visualization is exploratory and should not be interpreted as direct evidence that a specific Leiden cluster corresponds to a clinically defined GBM subtype.

---

## 12. Spatially Variable Gene Detection

SpatialDE is used to identify genes whose expression demonstrates statistically structured spatial variation.

SpatialDE addresses a different question from conventional differential-expression analysis:

```text
Differential Expression:
Which genes differ between transcriptomic clusters?

SpatialDE:
Which genes display non-random spatial expression patterns?
```

Because running SpatialDE across all retained genes is computationally expensive, the analysis is restricted to the top **2,000 highly variable genes**:

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000
)

adata_svg = adata[
    :,
    adata.var["highly_variable"]
].copy()
```

The resulting SpatialDE input contains approximately:

```text
2,870 spots × 2,000 genes
```

---

## 13. SpatialDE Input Preparation

The expression matrix is converted into a Pandas DataFrame:

```python
counts = pd.DataFrame(
    adata_svg.X.toarray(),
    columns=adata_svg.var_names,
    index=adata_svg.obs_names
)
```

The spatial coordinates are extracted using:

```python
coord = pd.DataFrame(
    adata_svg.obsm["spatial"],
    columns=[
        "x_coord",
        "y_coord"
    ],
    index=adata_svg.obs_names
)
```

Although the variable is named `counts`, the matrix contains normalized and log-transformed expression values rather than raw counts.

---

## 14. SpatialDE Analysis

SpatialDE is executed as:

```python
results = SpatialDE.run(
    coord,
    counts
)
```

SpatialDE is a statistical spatial-expression analysis method and is not a deep-learning model.

It does not use:

- CNNs,
- neural networks,
- transformers,
- or graph neural networks.

Its computational cost mainly comes from fitting and comparing spatial statistical models across many genes.

The results are integrated into the gene metadata:

```python
results = results.set_index(
    "g",
    drop=False
)

adata.var = adata.var.join(
    results.drop(columns="g"),
    how="left"
)
```

Candidate spatially variable genes are ranked by q-value:

```python
results.sort_values(
    "qval"
).head(10)
```

Lower q-values indicate stronger statistical evidence for spatially structured expression.

Selected genes such as:

```text
SLC9A3R2
COX20
```

are subsequently visualized on the tissue.

---

## 15. Literature-Guided GBM Gene Analysis

A previous study by Tang et al. used XGBoost to predict three major glioblastoma molecular subtypes:

- Classical,
- Mesenchymal,
- Proneural.

The study reported five genes with high predictive importance:

```text
NKAIN1
UBE2E2
F13A1
RNF149
PLAUR
```

In this project, these genes are **not used to train an XGBoost classifier**.

Instead, their spatial expression patterns are explored:

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color=[
        "NKAIN1",
        "UBE2E2",
        "F13A1",
        "RNF149",
        "PLAUR"
    ],
    alpha=0.6,
    cmap="plasma"
)
```

This provides a tissue-level view of genes previously associated with GBM subtype prediction.

---

## Main Methods

| Analysis Step | Method |
|---|---|
| Data structure | AnnData |
| Spatial technology | 10x Visium |
| Quality control | Scanpy |
| Normalization | Total-count normalization |
| Transformation | Log1p |
| Feature selection | Highly Variable Genes |
| Dimensionality reduction | PCA |
| Visualization | UMAP |
| Graph construction | Nearest-neighbor graph |
| Clustering | Leiden |
| Differential expression | t-test |
| Spatial visualization | Scanpy |
| Spatial gene detection | SpatialDE |
| Literature-guided analysis | GBM-associated genes |

---

## Main Findings

The analysis indicates that:

- glioblastoma tissue exhibits substantial transcriptomic heterogeneity,
- multiple transcriptomically distinct spatial spot clusters can be detected,
- several clusters form spatially coherent tissue regions,
- sequencing depth and detected-gene counts vary spatially,
- cluster-associated genes can be identified using differential-expression analysis,
- EGFR expression shows heterogeneous spatial distribution,
- SpatialDE identifies genes with non-random spatial expression patterns,
- and previously reported GBM-related genes also exhibit heterogeneous expression across the tissue.

---

## Important Interpretation Notes

### Visium Spots Are Not Single Cells

Each Visium spot can contain transcripts from multiple neighboring cells.

Therefore, the identified clusters represent groups of spatial spots rather than individual cells.

### Leiden Clusters Are Not Automatically Cell Types

Additional biological annotation, marker-gene analysis, reference mapping, or deconvolution would be required to assign reliable cell-type identities.

### Clusters Are Not GBM Molecular Subtypes

The discovered transcriptomic clusters should not automatically be labeled as classical, mesenchymal, or proneural GBM.

### SpatialDE Was Restricted to 2,000 Genes

SpatialDE was applied only to the top 2,000 highly variable genes rather than the full retained transcriptome.

This substantially reduces computational cost but may exclude genes with lower overall variability that still display spatial structure.

### XGBoost Was Not Trained in This Project

The XGBoost analysis belongs to previously published literature.

The five genes reported in that study are used here only for spatial visualization and exploratory biological interpretation.

---

## Computational Considerations

The most computationally demanding part of the project is:

```python
SpatialDE.run(coord, counts)
```

The runtime is caused by repeated spatial statistical model fitting across many genes rather than neural-network training.

Restricting the analysis from approximately:

```text
18,618 genes
```

to:

```text
2,000 highly variable genes
```

substantially reduces runtime and memory usage.

---

## Dependencies

The main Python packages used in this project are:

```text
scanpy
anndata
pandas
numpy
matplotlib
seaborn
igraph
leidenalg
tqdm
SpatialDE
```

For modern Python and SciPy environments, SpatialDE can be installed using:

```bash
pip install spatialde-modern
```

while retaining the standard import:

```python
import SpatialDE
```

---

## Installation

```bash
pip install scanpy
pip install pandas numpy matplotlib seaborn
pip install igraph leidenalg
pip install tqdm
pip install spatialde-modern
```

---

## Future Improvements

Potential extensions of this project include:

- biological annotation of the Leiden clusters,
- cell-type deconvolution,
- pathway and gene-set enrichment analysis,
- Moran's I and other spatial autocorrelation statistics,
- spatial neighborhood enrichment,
- ligand-receptor interaction analysis,
- GBM subtype signature scoring,
- analysis of multiple glioblastoma tissue sections,
- comparison across patients,
- and integration with histopathological image features.

---

## Project Scope

This repository should primarily be considered a:

> **Computational spatial transcriptomics analysis of glioblastoma tissue**

The project combines:

- transcriptomic quality control,
- unsupervised dimensionality reduction,
- graph-based clustering,
- spatial visualization,
- differential gene-expression analysis,
- spatially variable gene detection,
- and literature-guided biological exploration.

It is not currently a deep-learning or supervised subtype-classification project.

---

## Reference

Tang Y, Qazi MA, Brown KR, Mikolajewicz N, Moffat J, Singh SK, McNicholas PD.

**Identification of five important genes to predict glioblastoma subtypes.**

*Neuro-Oncology Advances.* 2021;3(1):vdab144.

DOI: `10.1093/noajnl/vdab144`

PMID: `34765972`

PMCID: `PMC8577514`

---

## Disclaimer

This project is intended for research and educational purposes only.

The analyses are exploratory and should not be interpreted as clinical diagnosis, validated tumor classification, or patient-specific medical guidance.

---

## Author

**Shayan Rokhva**

Research interests include:

- Artificial Intelligence
- Machine Learning
- Deep Learning
- Computational Biology
- Medical Image Analysis
- Spatial Transcriptomics / Omics
- Biomedical Data Science

---
shayanrokhva1999@gmail.com

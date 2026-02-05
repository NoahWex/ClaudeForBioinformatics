# Workshop Data

This directory uses a **manifest-based data reference system** - no data files are copied here.

## Data Manifest

The workshop data is defined in `workshop_data.yaml`, which references files at their source location.

**Why manifest-based?**
- No duplicate copies of large files
- Single source of truth for data locations
- Works on both HPC and local CRSP mount
- Consistent with lab conventions

## Data Source

**iHBCA Component Studies** from Reed et al. (2024), Nature Genetics.
- scVI-projected embeddings from the integrated Human Breast Cell Atlas
- Consistent annotation framework across studies

Reed AD, et al. "A human breast atlas integrating single-cell proteomics and transcriptomics." Nature Genetics, 2024.

## Files

| File | Size | Cells | Description |
|------|------|-------|-------------|
| `gray.rds` | ~255 MB | 52,681 | Gray et al. breast epithelial cells |
| `pal.rds` | ~1.8 GB | 80,000 | Pal et al. breast tissue atlas |

## Data Location

Data is accessed via manifest at:
```
HPC:   /share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/
Local: [CRSP mount]/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/
```

The `manifest_utils.R` functions automatically detect the environment and use the correct path.

## Data Provenance

### Gray et al.
- **Publication**: Gray GK, et al. "A human breast atlas..." (part of iHBCA)
- **Tissue**: Human breast
- **Technology**: 10x Genomics scRNA-seq
- **Cell types**: Epithelial subtypes (Luminal, Basal, Myoepithelial)
- **Key metadata columns**: `cell_type`, `sample_id`

### Pal et al.
- **Publication**: Pal B, et al. "A single-cell RNA expression atlas..." (part of iHBCA)
- **Tissue**: Human breast
- **Technology**: 10x Genomics scRNA-seq
- **Cell types**: Full breast tissue including immune and stromal
- **Key metadata columns**: `celltype_major`, `celltype_minor`, `SampleID`

## Checksums

For data integrity verification:

```bash
# Generate checksums
md5sum raw/*.rds > checksums.md5

# Verify checksums
md5sum -c checksums.md5
```

## Usage Notes

1. **Load in R**: `seurat_obj <- readRDS("raw/gray.rds")`
2. **Memory requirements**: ~4GB RAM for gray, ~16GB for pal
3. **Format**: Seurat v4/v5 objects with RNA assay

## Workshop Usage

- **Day 1**: Load and explore `gray.rds`
- **Day 2**: Compare metadata between gray and pal
- **Day 3**: Process both on HPC
- **Day 4**: Integrate datasets

## License

Data used under academic research provisions. Please cite the original publication.

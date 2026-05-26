---
description: "Debug an HDF5 file issue in metacooler. Inspects file structure, validates schema, checks for common problems like dtype mismatches or broken soft links."
---
# Debug an HDF5 file issue

Investigate the HDF5 issue described below. Follow these steps:

1. **Inspect the file** — use `uv run metacooler info <file>`, `uv run metacooler tree <file>`, `uv run metacooler attrs <file>` to understand the current state
2. **Check schema** — run `uv run python -c "from metacooler.schema import validate; validate('<file>')"`
3. **Common issues to check**:
   - String dtype (`object`/`StringDtype`) written directly to HDF5 instead of fixed-length bytes
   - Missing soft link `chroms → contigs`
   - Pixel ordering (bin1_id > bin2_id violates upper-triangle)
   - Index array length mismatches (`chrom_offset` vs n_contigs+1)
   - `genome_id` values out of range
4. **Report** the root cause and suggest a fix

## Issue: {{input}}

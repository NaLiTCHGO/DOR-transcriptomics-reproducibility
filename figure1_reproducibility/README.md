# DOR Figure 1 final reproducibility package

Public source package: repository release `v1.0.1`  
Status: `ACCEPTED GENERATOR / REAL REGENERATION PASS`

## One-command generation

Run from this package root:

```powershell
python code/generate_Figure1.py --package-root .
```

The project-specific generator reads `config/Figure1_layout_spec.json`, verifies the frozen provenance inputs, calls only the renderer and validators frozen under `code/vendor/`, and generates every final visual file under `outputs/`.

## Required software

- Python 3.12 or compatible
- Pillow 12.x
- Inkscape 1.4.x at `C:\Program Files\Inkscape\bin\inkscape.exe`, or pass another executable through `--inkscape`
- Arial and Arial Bold, or the documented Liberation/DejaVu fallbacks

## Final outputs

- `outputs/DOR_Figure1_final_v01.svg`: editable vector artwork with a 180-mm physical width
- `outputs/DOR_Figure1_final_v01.pdf`: one-page vector PDF with embedded Arial fonts
- `outputs/DOR_Figure1_final_v01_600dpi.png`: RGB 600-dpi PNG, 4252 × 4274 px
- `outputs/DOR_Figure1_final_v01_600dpi.tiff`: RGB LZW TIFF, 4252 × 4274 px at 600 dpi
- `outputs/DOR_Figure1_final_v01_preview.png`: review and manuscript-placement preview
- `outputs/Figure1_legend_final.txt`: final English legend

## Reproducibility chain

`frozen upstream summaries -> Figure1_layout_spec.json -> code/generate_Figure1.py -> frozen renderer -> SVG -> Inkscape PDF/PNG -> Pillow TIFF and DPI normalization`

The workflow is deterministic (`stochastic=false`). No manual post-edit is required. The source manifest records every frozen scientific source used to support counts and wording.

Before public release, the package was copied to an independent temporary directory and regenerated with the documented one-command workflow. Preview PNG, 600-dpi PNG, and LZW TIFF were pixel-identical to the accepted locked outputs; SVG and legend text were byte-identical.

The repository distribution intentionally omits generated `outputs/` and `qc/` files. After running the generator once, deterministic regeneration can be checked with:

```powershell
python code/run_reproducibility_audit.py --package-root .
```

The separate visual-equivalence validator additionally requires an accepted preview supplied with `--approved-preview`.

## Manuscript insertion

Insert at 180 mm when the journal permits a full-width two-column figure. Preserve the aspect ratio. Do not resize individual formats independently or edit text in a graphics program; revise the layout spec and rerun the generator instead.

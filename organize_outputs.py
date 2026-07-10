from pathlib import Path
import shutil
import math

PROJECT = Path.cwd()
PLOTS = PROJECT / "plots"
SORTED = PROJECT / "plots_sorted"

FOLDERS = {
    "01_exploratory_PCA_distance": [
        "PCA_Tumor_vs_Control.*",
        "sample_distance_heatmap.*"
    ],
    "02_volcano_and_DEG_heatmaps": [
        "volcano_plot.*",
        "volcano_gene_symbols.*",
        "top50_DEG_heatmap_ordered.*",
        "tumor_promoting_genes_heatmap.*",
        "tumor_opposing_genes_heatmap.*"
    ],
    "03_GO_enrichment": ["GO_*"],
    "04_KEGG_enrichment": ["KEGG_*"],
    "05_KEGG_pathview_maps": [],
    "06_Reactome_enrichment": ["Reactome_*"],
    "07_GSEA_Hallmark": ["GSEA_Hallmark_*"],
    "08_STRING_PPI_network": ["STRING_network*"],
    "09_hub_gene_heatmaps": ["Hub_Genes_Heatmap.*"],
    "10_biomarker_boxplots/individual_boxplots": [],
    "11_top_biomarker_heatmap": ["top_biomarker_heatmap.*"]
}

def copy_patterns():
    if not PLOTS.exists():
        print("[WARN] plots/ folder does not exist yet. Skipping output organization.")
        return

    if SORTED.exists():
        shutil.rmtree(SORTED)

    for folder, patterns in FOLDERS.items():
        out_dir = SORTED / folder
        out_dir.mkdir(parents=True, exist_ok=True)

        for pattern in patterns:
            for file in PLOTS.glob(pattern):
                if file.is_file():
                    shutil.copy2(file, out_dir / file.name)

    pathview_src = PLOTS / "pathview"
    pathview_dst = SORTED / "05_KEGG_pathview_maps"
    pathview_dst.mkdir(parents=True, exist_ok=True)
    if pathview_src.exists():
        for file in pathview_src.glob("*"):
            if file.is_file():
                shutil.copy2(file, pathview_dst / file.name)

    gsea_src = PLOTS / "GSEA_enrichment_plots"
    gsea_dst = SORTED / "07_GSEA_Hallmark" / "individual_GSEA_enrichment_plots"
    if gsea_src.exists():
        shutil.copytree(gsea_src, gsea_dst, dirs_exist_ok=True)

    boxplot_src = PLOTS / "hub_boxplots"
    boxplot_dst = SORTED / "10_biomarker_boxplots" / "individual_boxplots"
    boxplot_dst.mkdir(parents=True, exist_ok=True)
    if boxplot_src.exists():
        for file in boxplot_src.glob("*.png"):
            shutil.copy2(file, boxplot_dst / file.name)

def make_boxplot_collage():
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        print("[WARN] Pillow not installed. Skipping boxplot collage.")
        return

    boxplot_dir = SORTED / "10_biomarker_boxplots" / "individual_boxplots"
    output_file = SORTED / "10_biomarker_boxplots" / "all_biomarker_boxplots_landscape_collage.png"

    files = sorted(boxplot_dir.glob("*.png"))
    if not files:
        print("[WARN] No biomarker boxplot PNG files found. Skipping collage.")
        return

    thumb_w, thumb_h = 520, 420
    label_h = 60
    cols = 5
    rows = math.ceil(len(files) / cols)

    canvas = Image.new("RGB", (cols * thumb_w, rows * (thumb_h + label_h)), "white")
    draw = ImageDraw.Draw(canvas)

    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 22)
    except Exception:
        font = ImageFont.load_default()

    for i, file in enumerate(files):
        img = Image.open(file).convert("RGB")
        img.thumbnail((thumb_w - 30, thumb_h - 30))

        col = i % cols
        row = i // cols

        x = col * thumb_w + (thumb_w - img.width) // 2
        y = row * (thumb_h + label_h) + 10

        canvas.paste(img, (x, y))

        label = file.stem.replace("_boxplot", "")
        draw.text(
            (col * thumb_w + 20, row * (thumb_h + label_h) + thumb_h + 8),
            label,
            fill="black",
            font=font
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_file, quality=95)
    print(f"[OK] Biomarker boxplot collage created: {output_file}")

def organize_outputs():
    copy_patterns()
    make_boxplot_collage()
    print(f"[OK] Sorted plot folders created: {SORTED}")

if __name__ == "__main__":
    organize_outputs()

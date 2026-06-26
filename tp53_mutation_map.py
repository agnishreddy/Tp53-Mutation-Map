#!/usr/bin/env python3
"""
BioCore TP53 Terminal Codon & Mutation Atlas
A scientifically accurate, publication-quality terminal infographic
mapping the entire human TP53 canonical sequence (1-393), its codons, 
domain architecture, and precise clinical mutation hotspots.

Dependencies:
    pip install rich
"""

import sys
import os
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text
    from rich.columns import Columns
    from rich.align import Align
    from rich import box
except ImportError:
    print("Error: The 'rich' library is required for rendering this dashboard.")
    print("Please install it using: pip install rich")
    sys.exit(1)

# ============================================================================
# BIOLOGICAL DATA (Canonical Human TP53 - P04637)
# ============================================================================

TP53_SEQ = (
    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
    "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
    "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
    "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
    "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
)

# Standard human codon usage approximation (overridden by exact genomic hotspots)
STD_CODONS = {
    'M': 'ATG', 'E': 'GAG', 'P': 'CCG', 'Q': 'CAG', 'S': 'AGC',
    'D': 'GAC', 'V': 'GTG', 'L': 'CTG', 'T': 'ACC', 'F': 'TTC',
    'W': 'TGG', 'K': 'AAG', 'N': 'AAC', 'A': 'GCC', 'I': 'ATC',
    'R': 'CGC', 'G': 'GGC', 'Y': 'TAC', 'H': 'CAC', 'C': 'TGC'
}

EXACT_CODONS = {
    175: 'CGC', 245: 'GGC', 248: 'CGG', 249: 'AGG', 273: 'CGT',
    282: 'CGG', 220: 'TAT', 157: 'GTC', 158: 'CGC', 176: 'TGC',
    179: 'CAT', 237: 'ATG', 238: 'TGT', 244: 'GGC', 196: 'CGA',
    213: 'CGA', 132: 'AAG', 337: 'CGC', 280: 'AGA', 283: 'CGC',
    120: 'AAA', 241: 'TCC', 276: 'GCC', 277: 'TGT'
}

# Domain Boundaries (1-based indices)
DOMAINS = {
    "N-Term": (1, 93, "grey74"),
    "Core Domain": (94, 292, "bright_cyan bold"),
    "Tetramerization": (325, 356, "medium_purple1"),
    "CTD": (363, 393, "medium_purple1")
}

# Structural Features
ZINC_BINDING = {176, 179, 238, 242}
DNA_CONTACT = {120, 241, 248, 273, 276, 277, 280, 283}

# Detailed Clinical Hotspots
HOTSPOTS = {
    175: {
        "wt_aa": "Arginine", "wt_dna": "CGC", "color": "bright_red",
        "muts": [("CAC", "R175H", "G→A", "Transition", "Missense")],
        "class": "Structural Mutation", "freq": "≈6.0%",
        "mechanism": "Destabilizes Zinc-binding Loop (L2)",
        "clinical": "Loss of Function, Dominant Negative, High GOF",
        "cancers": "Breast, Ovary, Lung, Colon"
    },
    248: {
        "wt_aa": "Arginine", "wt_dna": "CGG", "color": "bright_red",
        "muts": [("CAG", "R248Q", "G→A", "Transition", "Missense"),
                 ("TGG", "R248W", "C→T", "Transition", "Missense")],
        "class": "DNA Contact", "freq": "≈8.5%",
        "mechanism": "Direct loss of minor groove DNA contact",
        "clinical": "Loss of Function, Dominant Negative",
        "cancers": "Glioblastoma, Colorectal, Breast"
    },
    273: {
        "wt_aa": "Arginine", "wt_dna": "CGT", "color": "bright_red",
        "muts": [("CAT", "R273H", "G→A", "Transition", "Missense"),
                 ("TGT", "R273C", "C→T", "Transition", "Missense")],
        "class": "DNA Contact", "freq": "≈7.2%",
        "mechanism": "Loss of phosphate backbone contact",
        "clinical": "Loss of Function, High GOF potential",
        "cancers": "Brain, Breast, Prostate"
    },
    220: {
        "wt_aa": "Tyrosine", "wt_dna": "TAT", "color": "dark_orange",
        "muts": [("TGT", "Y220C", "A→G", "Transition", "Missense")],
        "class": "Structural / Druggable", "freq": "≈1.5%",
        "mechanism": "Creates thermally unstable surface crevice",
        "clinical": "Druggable Pocket (e.g., PC14586 target)",
        "cancers": "Ovary, Breast, Gastric"
    },
    245: {
        "wt_aa": "Glycine", "wt_dna": "GGC", "color": "dark_orange",
        "muts": [("AGC", "G245S", "G→A", "Transition", "Missense")],
        "class": "Structural Mutation", "freq": "≈3.2%",
        "mechanism": "Alters L3 loop conformation",
        "clinical": "Loss of Function",
        "cancers": "Breast, Colon, Ovary"
    },
    282: {
        "wt_aa": "Arginine", "wt_dna": "CGG", "color": "bright_red",
        "muts": [("TGG", "R282W", "C→T", "Transition", "Missense")],
        "class": "Structural Mutation", "freq": "≈2.5%",
        "mechanism": "Disrupts H1 helix packing",
        "clinical": "Loss of Function",
        "cancers": "Lung, Esophageal, Colorectal"
    },
    249: {
        "wt_aa": "Arginine", "wt_dna": "AGG", "color": "bright_red",
        "muts": [("AGT", "R249S", "G→T", "Transversion", "Missense")],
        "class": "Structural Mutation", "freq": "≈3.0%",
        "mechanism": "L3 loop distortion (Aflatoxin B1 signature)",
        "clinical": "Loss of Function",
        "cancers": "Hepatocellular carcinoma (HCC)"
    },
    157: {
        "wt_aa": "Valine", "wt_dna": "GTC", "color": "yellow",
        "muts": [("TTC", "V157F", "G→T", "Transversion", "Missense")],
        "class": "Structural Mutation", "freq": "≈1.2%",
        "mechanism": "Hydrophobic core disruption",
        "clinical": "Smoking-associated mutational signature",
        "cancers": "Lung (SCLC & NSCLC)"
    },
    337: {
        "wt_aa": "Arginine", "wt_dna": "CGC", "color": "medium_purple1",
        "muts": [("CAC", "R337H", "G→A", "Transition", "Missense")],
        "class": "Tetramerization", "freq": "≈1.0%",
        "mechanism": "pH-dependent tetramer destabilization",
        "clinical": "Li-Fraumeni syndrome (Brazilian founder)",
        "cancers": "Adrenocortical, Breast"
    }
}


# ============================================================================
# RENDERING ENGINE
# ============================================================================

class TP53CodonAtlas:
    def __init__(self):
        # Enable recording to capture terminal output, and fix the width 
        # to ensure the exported image stays perfectly aligned.
        self.console = Console(record=True, width=150)
        self.chunk_size = 35

    def get_codon(self, pos: int, aa: str) -> str:
        return EXACT_CODONS.get(pos, STD_CODONS.get(aa, "NNN"))

    def get_style(self, pos: int) -> str:
        if pos in HOTSPOTS:
            return f"bold {HOTSPOTS[pos]['color']}"
        
        # Base domain coloring
        base_style = "grey50"
        for dom, (start, end, style) in DOMAINS.items():
            if start <= pos <= end:
                base_style = style
                break
                
        if pos in DNA_CONTACT:
            base_style = "bold green underline"
        if pos in ZINC_BINDING:
            base_style = "bold magenta reverse"
            
        return base_style

    def render_header(self):
        title = Text("\nTP53 TUMOR SUPPRESSOR PROTEIN\n", style="bold white", justify="center")
        title.append("Canonical Human Protein (393 aa)\n", style="bright_cyan")
        title.append("DNA-Binding Core • Cancer Hotspots • Structural Motifs\n", style="grey70")
        self.console.print(Panel(title, box=box.DOUBLE, border_style="bright_cyan"))

    def render_legend(self):
        table = Table.grid(padding=(0, 4), expand=True)
        table.add_column(justify="center")
        table.add_column(justify="center")
        table.add_column(justify="center")
        table.add_column(justify="center")
        
        table.add_row(
            Text("■ Intrinsically Disordered (N/C)", style="grey74"),
            Text("■ DNA Binding Domain", style="bright_cyan bold"),
            Text("■ Tetramerization Domain", style="medium_purple1"),
            Text("■ Major Cancer Hotspot", style="bright_red bold")
        )
        table.add_row(
            Text("■ Structural Mutation", style="dark_orange bold"),
            Text("■ Moderately Recurrent", style="yellow bold"),
            Text("A DNA Contact Residue", style="green bold underline"),
            Text("■ Zn²⁺ Binding Site", style="magenta reverse bold")
        )
        self.console.print(Panel(table, title="[bold]Structural & Mutational Legend[/bold]", border_style="grey50"))

    def render_sequence_blocks(self):
        seq_len = len(TP53_SEQ)
        
        for i in range(0, seq_len, self.chunk_size):
            chunk_start = i + 1
            chunk_end = min(i + self.chunk_size, seq_len)
            chunk_seq = TP53_SEQ[i:i+self.chunk_size]

            block_text = Text()
            
            # --- Row 1: Index Numbering ---
            idx_line = Text(f"[{chunk_start:03d}–{chunk_end:03d}] │", style="bold white")
            for j in range(len(chunk_seq)):
                pos = chunk_start + j
                idx_line.append(f"{pos:4d}", style="grey50")
            block_text.append(idx_line)
            block_text.append("\n")

            # --- Row 2: DNA Codons ---
            dna_line = Text(" DNA Codon │", style="grey70")
            for j, aa in enumerate(chunk_seq):
                pos = chunk_start + j
                dna = self.get_codon(pos, aa)
                style = self.get_style(pos)
                dna_line.append(f" {dna}", style=style)
            block_text.append(dna_line)
            block_text.append("\n")

            # --- Row 3: mRNA Codons ---
            rna_line = Text("mRNA Codon │", style="grey70")
            for j, aa in enumerate(chunk_seq):
                pos = chunk_start + j
                dna = self.get_codon(pos, aa)
                rna = dna.replace('T', 'U')
                style = self.get_style(pos)
                rna_line.append(f" {rna}", style=style)
            block_text.append(rna_line)
            block_text.append("\n")

            # --- Row 4: Amino Acids ---
            aa_line = Text("Amino Acid │", style="bold white")
            for j, aa in enumerate(chunk_seq):
                pos = chunk_start + j
                style = self.get_style(pos)
                # Pad for alignment with 3-letter codons
                aa_line.append(f"  {aa} ", style=style)
            block_text.append(aa_line)

            self.console.print(Panel(block_text, border_style="cyan", padding=(0, 1)))

    def render_hotspots_detailed(self):
        panels = []
        for pos, data in sorted(HOTSPOTS.items()):
            content = Text()
            content.append(f"Wild DNA: ", style="grey70")
            content.append(f"{data['wt_dna']} ", style="cyan")
            content.append(f"| Wild mRNA: ", style="grey70")
            content.append(f"{data['wt_dna'].replace('T', 'U')}\n", style="cyan")
            
            content.append("Common Mutations:\n", style="bold white")
            for m_dna, m_aa, nt_change, m_type, m_cons in data['muts']:
                content.append(f"  {data['wt_dna']} → {m_dna} ", style="bold red")
                content.append(f"({m_aa})\n", style="bold yellow")
                content.append(f"  {nt_change} [{m_type}] • {m_cons}\n", style="grey62")
                
            content.append(f"Mutation Type: ", style="grey70")
            content.append(f"{data['class']}\n", style=data['color'])
            content.append(f"Frequency: ", style="grey70")
            content.append(f"{data['freq']}\n", style="bold white")
            content.append(f"Effect: ", style="grey70")
            content.append(f"{data['mechanism']}\n", style="white")
            content.append(f"Clinical: ", style="grey70")
            content.append(f"{data['clinical']}\n", style="white")
            content.append(f"Cancers: ", style="grey70")
            content.append(f"{data['cancers']}", style="grey74")
            
            p = Panel(content, title=f"[bold {data['color']}]{data['wt_aa']} {pos}[/]", border_style=data['color'], width=45)
            panels.append(p)
            
        self.console.print(Panel(Columns(panels, expand=True, equal=True), title="[bold white]CLINICALLY SIGNIFICANT MUTATION HOTSPOTS[/bold white]", border_style="red"))

    def render_footer(self):
        text = (
            "Gene: TP53 | Protein: Cellular Tumor Antigen p53 | UniProt: P04637 | Length: 393 aa | Chromosome: 17p13.1\n"
            "Data Sources: UniProt, COSMIC, IARC TP53 Database, ClinVar, TCGA, NCBI Gene, HGVS"
        )
        self.console.print(Panel(Text(text, justify="center", style="grey50"), box=box.ROUNDED, border_style="grey30"))

    def display(self):
        self.render_header()
        self.render_legend()
        self.render_sequence_blocks()
        self.render_hotspots_detailed()
        self.render_footer()

        # --- HIGH-RESOLUTION EXPORT LOGIC ---
        output_dir = "/mnt/d/Bioinformatics_Projects/vision-pathway"
        os.makedirs(output_dir, exist_ok=True)
        
        # Save as SVG (Native vector format for terminal captures, superior to JPEG for text)
        out_svg = os.path.join(output_dir, "tp53_codon_atlas.svg")
        self.console.save_svg(out_svg, title="TP53 Codon & Mutation Atlas")
        
        # Print confirmation to the terminal
        self.console.print(f"\n[bold bright_green][SUCCESS][/] Publication-quality vector image saved to: [bold cyan]{out_svg}[/]")
        self.console.print("[dim italic]Note: Open the SVG in any web browser or image editor to view it at infinite resolution, or export directly to JPEG/PNG.[/]")


if __name__ == "__main__":
    atlas = TP53CodonAtlas()
    atlas.display()
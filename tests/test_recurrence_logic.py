import subprocess
from pathlib import Path
import pandas as pd


def test_relaxed_vs_strict(tmp_path):
    rows = [
        ["G1","BRCA1",1,0.4,2], ["G1","BRCA2",1,0.5,2],
        ["G2","BRCA1",2,0.4,3], ["G2","BRCA2",2,0.6,3],
    ]
    inp = tmp_path / "in.tsv"
    pd.DataFrame(rows, columns=["GENE","GWAS","N_CLUMPS","MAX_ABS_BETA","N_MAPPED_SNPS"]).to_csv(inp, sep="\t", index=False)
    script = Path(__file__).resolve().parents[1] / "scripts" / "07_recurrence_summary.py"
    relaxed = tmp_path / "relaxed.tsv"
    strict = tmp_path / "strict.tsv"
    subprocess.run(["python", str(script), "--inputs", str(inp), "--min-gwas", "2", "--min-clumps", "1", "--output", str(relaxed)], check=True)
    subprocess.run(["python", str(script), "--inputs", str(inp), "--min-gwas", "2", "--min-clumps", "2", "--output", str(strict)], check=True)
    assert set(pd.read_csv(relaxed, sep="\t").GENE) == {"G1", "G2"}
    assert set(pd.read_csv(strict, sep="\t").GENE) == {"G2"}

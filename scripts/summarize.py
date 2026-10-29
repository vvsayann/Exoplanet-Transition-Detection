"""Summarize current experiment checkpoints without repeating archive requests."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import json
import pandas as pd
from src.catalog import compare_known


def main():
    for stage in ['known','blind','labeled']:
        directory = Path('results/research')/stage
        path = directory/'frozen_results.json'
        rows = json.loads(path.read_text()) if path.exists() else []
        errors = directory/'failures.json'
        failures = json.loads(errors.read_text()) if errors.exists() else []
        print(f'{stage}: {len(rows)} completed, {len(failures)} failed')
        if stage == 'known' and rows and Path('data/planets.csv').exists():
            compare_known(rows, pd.read_csv('data/planets.csv')).to_csv('results/research/known_comparison.csv',index=False)
        for row in rows:
            fit = row.get('fit',{}).get('diagnostics',{})
            print(f"  {row['target']}: P={row['period']:.8f}, SNR={row['depth_snr']:.1f}, converged={fit.get('converged', 'not fitted')}")


if __name__ == '__main__':
    main()

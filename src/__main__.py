import argparse
import json
from pathlib import Path
import numpy as np
import lightkurve as lk
from .pipeline import run_pipeline
from .experiments import synthetic_curve, injection_recovery, batch, evaluate_labeled
from .catalog import catalog, compare_known


def main():
    parser = argparse.ArgumentParser(description='Exoplanet transit detection research pipeline')
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ['run', 'demo', 'batch']:
        p = sub.add_parser(command)
        if command == 'run':
            p.add_argument('target'); p.add_argument('--mission', choices=['Kepler','TESS'], default='Kepler')
            p.add_argument('--stellar-radius', type=float); p.add_argument('--stellar-radius-error', type=float, default=0.)
        if command == 'batch':
            p.add_argument('manifest')
        p.add_argument('--output', default=f'results/{command}')
        p.add_argument('--no-fit', action='store_true')
        p.add_argument('--steps', type=int, default=1500); p.add_argument('--burn', type=int, default=500)
        p.add_argument('--max-period', type=float, default=50.); p.add_argument('--max-segments', type=int, default=4)
    p = sub.add_parser('catalog'); p.add_argument('table', choices=['pscomppars','cumulative','toi']); p.add_argument('--output', default='data/catalog.csv')
    p = sub.add_parser('compare'); p.add_argument('results'); p.add_argument('catalog'); p.add_argument('--output', default='results/comparison.csv')
    p = sub.add_parser('evaluate'); p.add_argument('results'); p.add_argument('labels'); p.add_argument('--output', default='results/vetting_evaluation.json')
    p = sub.add_parser('inject'); p.add_argument('lightcurve', help='Pipeline lightcurve.npz, or synthetic for an offline demonstration')
    p.add_argument('--periods', type=float, nargs='+', default=[1.,3.,7.]); p.add_argument('--radii', type=float, nargs='+', default=[2.,4.,8.])
    p.add_argument('--repeats', type=int, default=5); p.add_argument('--output', default='results/injection')
    args = vars(parser.parse_args()); command = args.pop('command')
    if command in ['run','demo','batch']:
        args['fit'] = not args.pop('no_fit')
        if command == 'demo':
            args['max_period'] = min(args['max_period'], 10.)
            result = run_pipeline('SYNTHETIC_DEMO', mission='SYNTHETIC', curves=[synthetic_curve()], stellar_radius=1., **args)
        elif command == 'run':
            result = run_pipeline(args.pop('target'), **args)
        else:
            result = batch(args.pop('manifest'), **args)
        print(json.dumps(result, indent=2))
    elif command == 'catalog':
        catalog(**args)
    elif command == 'compare':
        import pandas as pd
        rows = json.loads(Path(args['results']).read_text())
        compare_known(rows, pd.read_csv(args['catalog'])).to_csv(args['output'], index=False)
    elif command == 'evaluate':
        print(evaluate_labeled(args['results'], args['labels'], args['output']))
    elif command == 'inject':
        path = args.pop('lightcurve')
        if path == 'synthetic':
            # Quiet synthetic baseline, explicitly distinct from real-star completeness.
            rng = np.random.default_rng(42); t = np.arange(0,27,1/48)
            curves = [lk.LightCurve(time=t, flux=1+rng.normal(0,.0003,len(t)), flux_err=np.full(len(t),.0003))]
        else:
            data = np.load(path)
            curves = [lk.LightCurve(time=data['time'], flux=data['flux'], flux_err=data['flux_err'])]
        injection_recovery(curves, **args)


if __name__ == '__main__':
    main()

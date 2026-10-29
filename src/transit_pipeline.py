"""Compatibility entry point. Prefer `python -m src run Kepler-10`."""
if __name__ == '__main__':
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.__main__ import main
    if len(sys.argv) == 1:
        sys.argv.extend(['run', 'Kepler-10'])
    main()

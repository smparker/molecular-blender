"""Render ORCA ethene orbital 7 (Molden orbital 8).

blender -b --factory-startup --python-exit-code 1 --python examples/render_orca_ethene.py -- OUTPUT_DIR [ethene.molden.input]
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_ethene import render

args = sys.argv[sys.argv.index('--') + 1:]
render('orca', args[0], Path(args[1]).resolve() if len(args) > 1 else None)

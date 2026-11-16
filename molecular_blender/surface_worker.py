"""Private entry point run by a background Blender process, never a thread."""
from pathlib import Path
import json
import os
import pickle
import sys
import traceback


def main(directory):
    # Prefer the same addon installation/check-out as the launching process.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from molecular_blender.molecule import Molecule
    from molecular_blender.plotter import compute_surface_vertex_sets

    class Progress(dict):
        def __setitem__(self, key, value):
            super().__setitem__(key, value)
            temporary = directory / 'progress.tmp'
            temporary.write_text(json.dumps(self))
            os.replace(temporary, directory / 'progress.json')

    try:
        with (directory / 'request.pkl').open('rb') as stream:
            data, options = pickle.load(stream)
        molecule = Molecule.from_dict('surface-worker', data)
        result = compute_surface_vertex_sets(molecule, options, Progress())
        with (directory / 'result.pkl').open('wb') as stream:
            pickle.dump(result, stream, protocol=4)
    except Exception as exc:
        (directory / 'error.txt').write_text(str(exc))
        traceback.print_exc()
        raise


if __name__ == '__main__':
    main(Path(sys.argv[sys.argv.index('--') + 1]))

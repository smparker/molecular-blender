"""Private entry point run by a background Blender process, never a thread."""
from pathlib import Path
import json
import os
import pickle
import sys
import time
import traceback


def main(directory):
    # Prefer the same addon installation/check-out as the launching process.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from molecular_blender.molecular_data import Molecule
    from molecular_blender.surface_calculation import compute_surface_vertex_sets, iter_surface_frames
    from molecular_blender.surface_jobs import SavedSurfaceSource

    class Progress(dict):
        def __setitem__(self, key, value):
            super().__setitem__(key, value)
            temporary = directory / 'progress.tmp'
            temporary.write_text(json.dumps(self))
            os.replace(temporary, directory / 'progress.json')

    try:
        with (directory / 'request.pkl').open('rb') as stream:
            data, options, stream_frames = pickle.load(stream)
        if isinstance(data, SavedSurfaceSource):
            molecule = data.load(options)
        else:
            molecule = Molecule.from_dict('surface-worker', data)
        if stream_frames:
            frames = iter_surface_frames(molecule, options, Progress(), indexed=True)
            count = 0
            for batch in frames:
                # At most two ready files, even if Blender is busy rendering or drawing.
                while count >= 2 and (directory / f'frame-{count-2}.pkl').exists():
                    time.sleep(.025)
                temporary = directory / 'frame.tmp'
                with temporary.open('wb') as stream:
                    pickle.dump(batch, stream, protocol=4)
                temporary.replace(directory / f'frame-{count}.pkl')
                count += 1
                del batch
            (directory / 'frame-count').write_text(str(count))
        else:
            result = compute_surface_vertex_sets(molecule, options, Progress(), indexed=True)
            with (directory / 'result.pkl').open('wb') as stream:
                pickle.dump(result, stream, protocol=4)
    except Exception as exc:
        (directory / 'error.txt').write_text(str(exc))
        traceback.print_exc()
        raise


if __name__ == '__main__':
    main(Path(sys.argv[sys.argv.index('--') + 1]))

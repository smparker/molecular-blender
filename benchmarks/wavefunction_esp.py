"""Time preview, final and cached ESP on a displayed 0.001 e/bohr³ density surface.

blender -b --factory-startup --python benchmarks/wavefunction_esp.py -- file.molden [net-charge]
Run one molecule per process. Peak RSS includes Blender, surface generation and ESP.
In Blender 4.3.2 on arm64 macOS, water had 2,664 vertices (0.04 s preview,
0.04 s final, 216 MiB peak). Standard-normalized PBE0 phenanthroline had
13,890 vertices (0.55 s preview, 4.10 s final, 278 MiB peak). Final uses
the contracted density prepared during preview. Times exclude worker startup.
"""
from pathlib import Path
from resource import getrusage, RUSAGE_SELF
import sys
from time import perf_counter

import bpy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from molecular_blender import scene_panel
from molecular_blender.plotter import BlendMoleculeSetup, draw_surfaces
from molecular_blender.scene_state import read_record
from molecular_blender.scene_surfaces import surface_groups


def main(path, charge=0):
    scene_panel.register()
    model, styler, opts = BlendMoleculeSetup(bpy.context, str(Path(path).resolve()), volume='density',
        isovalues='.001', cumulative=False, resolution=.16, remesh=0, molecular_charge=charge)
    draw_surfaces(model, styler, bpy.context, opts)
    owner = model.scene_owner
    key, group = next(iter(surface_groups(owner, read_record(owner)).items()))
    bpy.context.view_layer.objects.active = owner
    obj, = group['objects']
    samples = []
    for label, mode in [('preview', 'preview'), ('final', 'final'), ('cached', 'final')]:
        start = perf_counter()
        assert bpy.ops.molecular_blender.compute_esp(group=key, quality=mode) == {'FINISHED'}
        elapsed = perf_counter()-start
        values = np.array([a.value for a in obj.data.attributes['mb_esp'].data])
        samples.append(values)
        rss = getrusage(RUSAGE_SELF).ru_maxrss / (1024**(2 if sys.platform == 'darwin' else 1))
        print(f'{label}: {len(values)} vertices, {elapsed:.3f} s, peak RSS {rss:.1f} MiB', flush=True)
    print(f'Preview RMS difference: {np.sqrt(np.mean((samples[0]-samples[1])**2)):.6f} Hartree/e')


if __name__ == '__main__':
    args = sys.argv[sys.argv.index('--')+1:]
    main(args[0], int(args[1]) if len(args)>1 else 0)

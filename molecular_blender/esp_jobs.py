"""Sample validated Molden ESP on displayed vertices, inside a surface worker."""
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import os
from pathlib import Path
import pickle
import tempfile

import numpy as np

from .surface_jobs import SavedSurfaceSource
from .wavefunction_esp import HARTREE_IN_VOLTS, WavefunctionESP


@lru_cache(maxsize=1)
def session_cache():
    return tempfile.TemporaryDirectory(prefix='molecular-esp-')


def _read_cached(path, key):
    try:
        with path.open('rb') as stream:
            if pickle.load(stream) == key:
                return pickle.load(stream)
    except FileNotFoundError:
        pass
    return None


def _write_cached(path, key, value):
    # All files are private to the Blender session. Publish complete entries atomically.
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        pending = Path(stream.name)
        try:
            pickle.dump(key, stream)
            pickle.dump(value, stream, protocol=4)
            stream.close()
            os.replace(pending, path)
        finally:
            pending.unlink(missing_ok=True)


@dataclass
class ESPRequest:
    source: SavedSurfaceSource
    convention: str
    charge: int
    quality: str
    units: str
    points: np.ndarray
    cache: Path

    @classmethod
    def from_record(cls, record, points, convention, charge, quality, units):
        return cls(SavedSurfaceSource.from_record(record['source'], record), convention, charge,
                   quality, units, points, Path(session_cache().name))

    def calculate(self, progress):
        if self.quality not in ('preview', 'final') or self.units not in ('Hartree/e', 'V'):
            raise ValueError('Choose an ESP quality and potential unit')
        digest = hashlib.sha256()
        with open(self.source.path, 'rb') as stream:
            for block in iter(lambda: stream.read(1024*1024), b''):
                digest.update(block)
        digest.update(pickle.dumps((self.source, self.charge, self.convention), protocol=4))
        model_key = digest.digest()
        self.points = np.ascontiguousarray(self.points, dtype=np.float64)
        if self.points.ndim != 2 or self.points.shape[1] != 3 or not np.isfinite(self.points).all():
            raise ValueError('ESP needs finite (N, 3) point coordinates')
        digest.update(self.points.view(np.uint8))
        value_key = digest.digest()
        value_path = self.cache / (self.quality + '.pkl')
        values = _read_cached(value_path, value_key)
        if values is None:
            model_path = self.cache / 'density.pkl'
            field = _read_cached(model_path, model_key)
            if field is None:
                data, molecule = self.source.read(dict(molden_convention=self.convention))
                field = WavefunctionESP(data, self.charge, molecule.orbitals)
                _write_cached(model_path, model_key, field)
            if self.quality == 'preview':
                # Sample an actual displayed vertex per 0.45 Å cell. Adjacent vertex
                # colors interpolate across the faces; no potential grid covers nuclei.
                cells = np.floor(self.points/.45).astype(np.int64)
                _, indices, inverse = np.unique(cells, axis=0, return_index=True, return_inverse=True)
                points = self.points[indices]
            else:
                points, inverse = np.unique(self.points, axis=0, return_inverse=True)
            result = np.empty(len(points))
            progress['total'] = len(points)
            for start in range(0, len(points), 64):
                end = min(start+64, len(points))
                result[start:end] = field.batch(points[start:end])
                progress['completed'] = end
            values = result[inverse.reshape(-1)]
            _write_cached(value_path, value_key, values)
        return values * (HARTREE_IN_VOLTS if self.units == 'V' else 1.)


@dataclass
class ContourESPRequests:
    requests: dict

    def calculate(self, progress):
        return {level: request.calculate(progress) for level, request in self.requests.items()}

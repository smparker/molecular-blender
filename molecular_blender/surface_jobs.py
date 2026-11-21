"""Isolated, cancellable surface jobs. No threads run inside interactive Blender.

The private temporary directory holds trusted parent/child pickle messages. The
child uses the same Blender executable (and thus Python/NumPy ABI) as the parent.
Terminating the process interrupts even a compiled numerical loop.
"""
import atexit
from dataclasses import dataclass
import json
from pathlib import Path
import pickle
import subprocess
import tempfile


_jobs = set()


@dataclass
class SavedSurfaceSource:
    """A sidebar job can reopen the source directly inside the worker."""
    path: str
    atoms: list
    orbitals: list
    datasets: list = None

    @classmethod
    def from_record(cls, path, record):
        atoms = [(a['index'], a['element'], tuple(a['position'])) for a in record['atoms']]
        fields = record.get('fields')
        return cls(path, atoms, record['orbitals'], fields['datasets'] if fields else None)

    def load(self, options):
        from .importers import molecule_from_file
        from .molecule import Molecule
        from .scene_state import orbital_metadata

        data = molecule_from_file(self.path, options)
        model = Molecule.from_dict('surface-worker', data)
        actual = [(a.index, a.element.symbol, tuple(a.position)) for a in model.atoms]
        if actual != self.atoms:
            raise ValueError('Source atom identities or coordinates have changed. Reimport the source file.')
        if self.datasets is not None:
            if model.volume is None or model.volume.labels != self.datasets or model.volume_trajectory:
                raise ValueError('Source cube datasets or frame count have changed. Reimport the source file.')
        elif orbital_metadata(model.orbitals) != self.orbitals:
            raise ValueError('Source orbital metadata has changed. Reimport the Molden file.')
        return model


class SurfaceJob:
    def __init__(self, executable, molecule_data, options, *, stream_frames=False, on_close=None):
        self._temporary = tempfile.TemporaryDirectory(prefix='molecular-surface-')
        self.directory = Path(self._temporary.name)
        self._next_frame = 0
        self._on_close = on_close
        self.stream_frames = stream_frames
        self.process = None
        self.cancelled = False
        self._log = None
        try:
            with (self.directory / 'request.pkl').open('wb') as stream:
                pickle.dump((molecule_data, options, stream_frames), stream, protocol=4)
            self._log = (self.directory / 'worker.log').open('wb')
            self.process = subprocess.Popen([
                executable, '--background', '--factory-startup', '--disable-autoexec',
                '--threads', '1', '--python-exit-code', '1', '--python',
                str(Path(__file__).with_name('surface_worker.py')),
                '--', str(self.directory),
            ], stdout=self._log, stderr=subprocess.STDOUT)
            _jobs.add(self)
        except BaseException:
            self.close()
            raise

    @property
    def running(self):
        return self.process is not None and self.process.poll() is None

    def progress(self):
        try:
            return json.loads((self.directory / 'progress.json').read_text())
        except (OSError, ValueError):
            return {'completed': 0, 'total': 1}

    def take_frame(self):
        """Consume one atomically published batch and free its slot in the worker queue."""
        path = self.directory / f'frame-{self._next_frame}.pkl'
        try:
            stream = path.open('rb')
        except FileNotFoundError:
            return None
        with stream:
            batch = pickle.load(stream)
        path.unlink()
        self._next_frame += 1
        return batch

    def check_finished(self):
        if self.cancelled:
            raise RuntimeError('Surface calculation was cancelled')
        if self.running:
            raise RuntimeError('Surface calculation is still running')
        error_path = self.directory / 'error.txt'
        if error_path.exists():
            raise RuntimeError(error_path.read_text())
        if self.process is None or self.process.returncode != 0:
            log = (self.directory / 'worker.log').read_text(errors='replace')
            raise RuntimeError('Surface worker failed: ' + log[-2000:])
        if self.stream_frames and self._next_frame != int((self.directory / 'frame-count').read_text()):
            raise RuntimeError('Incomplete trajectory surface delivery')

    def result(self):
        self.check_finished()
        with (self.directory / 'result.pkl').open('rb') as stream:
            return pickle.load(stream)

    def close(self):
        """Reap the worker before removing its files; safe on every exit path."""
        if self.running:
            self.cancelled = True
            self.process.terminate()
            try:
                self.process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self._log is not None:
            self._log.close()
            self._log = None
        self._temporary.cleanup()
        _jobs.discard(self)
        if self._on_close is not None:
            callback, self._on_close = self._on_close, None
            callback()

    def cancel(self):
        self.cancelled = True
        self.close()


def cancel_all_jobs(*_args):
    for job in tuple(_jobs):
        job.cancel()


atexit.register(cancel_all_jobs)

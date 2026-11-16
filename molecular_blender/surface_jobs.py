"""Isolated, cancellable surface jobs. No threads run inside interactive Blender.

The private temporary directory holds trusted parent/child pickle messages. The
child uses the same Blender executable (and thus Python/NumPy ABI) as the parent.
Terminating the process interrupts even a compiled numerical loop.
"""
import atexit
import json
from pathlib import Path
import pickle
import subprocess
import tempfile


_jobs = set()


class SurfaceJob:
    def __init__(self, executable, molecule_data, options):
        self._temporary = tempfile.TemporaryDirectory(prefix='molecular-surface-')
        self.directory = Path(self._temporary.name)
        self.process = None
        self.cancelled = False
        self._log = None
        try:
            with (self.directory / 'request.pkl').open('wb') as stream:
                pickle.dump((molecule_data, options), stream, protocol=4)
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

    def result(self):
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

    def cancel(self):
        self.cancelled = True
        self.close()


def cancel_all_jobs(*_args):
    for job in tuple(_jobs):
        job.cancel()


atexit.register(cancel_all_jobs)

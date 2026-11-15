"""Installer checks use temporary addon trees and never touch Blender settings."""

from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import blender_setup


class TestInstallTargets(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / 'source'
        (self.source / 'molecular_blender').mkdir(parents=True)
        (self.source / 'molecular_blender' / 'version').write_text('new')
        self.versions = [self.root / version for version in ('4.3', '4.10')]
        for path in self.versions:
            addon = path / 'scripts/addons/molecular_blender'
            addon.mkdir(parents=True)
            (addon / 'version').write_text('old')
        self.addCleanup(patch.stopall)
        patch.object(blender_setup, 'currentDirectory', str(self.source)).start()
        patch.object(blender_setup, 'find_user_blender_dirs',
                     return_value=[str(path) for path in self.versions]).start()
        self.output = io.StringIO()

    def test_selected_version_only(self):
        with redirect_stdout(self.output):
            blender_setup.install(link_target='4.10')
        self.assertEqual((self.versions[0] / 'scripts/addons/molecular_blender/version').read_text(), 'old')
        self.assertEqual((self.versions[1] / 'scripts/addons/molecular_blender/version').read_text(), 'new')
        self.assertIn('replaces existing installation', self.output.getvalue())

    def test_dry_run_and_invalid_targets_preserve_installations(self):
        with redirect_stdout(self.output):
            targets = blender_setup.install(link_target='4.3,4.10', dry_run=True)
        self.assertEqual(len(targets), 2)
        with self.assertRaisesRegex(ValueError, 'Unknown Blender install target'):
            blender_setup.install(link_target='4.3,99.0')
        for path in self.versions:
            self.assertEqual((path / 'scripts/addons/molecular_blender/version').read_text(), 'old')

    @unittest.skipIf(os.name == 'nt', 'symlink creation requires Windows privileges')
    def test_replace_dangling_link(self):
        target = self.versions[0] / 'scripts/addons/molecular_blender'
        import shutil
        shutil.rmtree(target)
        target.symlink_to(self.root / 'missing', target_is_directory=True)
        with redirect_stdout(self.output):
            blender_setup.install(install_mode='link', link_target='4.3')
        self.assertTrue(target.is_symlink())
        self.assertEqual((target / 'version').read_text(), 'new')

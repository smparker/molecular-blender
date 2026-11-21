"""Assemble a Blender add-on from the same native binaries as the release wheels."""
from pathlib import Path
import shutil
import subprocess
import tempfile
from zipfile import ZipFile

from .addon_files import iterRelativeAddonFiles
from .export import execute_Export
from .generic import getAddonVersion


def export_wheels(source, destination, wheel_dir):
    source = Path(source)
    wheels = sorted(Path(wheel_dir).glob('molecular_blender-*.whl'))
    if not wheels:
        raise ValueError(f'No molecular_blender wheels in {wheel_dir}')
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        package = root / source.name
        output = Path(destination).resolve()
        for name in iterRelativeAddonFiles(source):
            src = source / name
            if src.suffix in {'.so', '.pyd', '.dll', '.dylib'} or src.resolve() == output:
                continue
            target = package / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, target)

        modules = {path.stem for path in source.glob('*.pyx')}
        version = '.'.join(map(str, getAddonVersion(source / '__init__.py')))
        inventory = {}
        mac = {}
        for wheel in wheels:
            _, release, python, abi, platform = wheel.stem.split('-')
            if platform.startswith('musllinux'):
                # Official Blender Linux downloads use glibc.
                continue
            if release != version or python != abi:
                raise ValueError(f'Expected version {version} and a CPython-specific ABI: {wheel.name}')
            if platform.startswith('macosx_') and platform.endswith(('_arm64', '_x86_64')):
                family = 'mac-arm64' if platform.endswith('_arm64') else 'mac-x86_64'
            elif platform.startswith('manylinux') and platform.endswith('_x86_64'):
                family = 'linux-x64'
            elif platform == 'win_amd64':
                family = 'windows-x64'
            else:
                raise ValueError(f'Unsupported add-on wheel platform: {wheel.name}')
            built = inventory.setdefault(family, set())
            if abi in built:
                raise ValueError(f'Duplicate {family} wheel for {abi}')
            built.add(abi)
            with ZipFile(wheel) as archive:
                found = {Path(n).name.split('.')[0] for n in archive.namelist()
                         if Path(n).parent == Path(source.name) and n.endswith(('.so', '.pyd'))}
                if found != modules:
                    raise ValueError(f'Expected {sorted(modules)} in {wheel.name}, found {sorted(found)}')
                if family == 'linux-x64':
                    _extract_linux(archive, package)
                elif family == 'windows-x64':
                    # The Windows extension loader searches beside the .pyd for its DLLs.
                    for name in archive.namelist():
                        if name.endswith(('.pyd', '.dll')):
                            (package / Path(name).name).write_bytes(archive.read(name))
                else:
                    arch = family.removeprefix('mac-')
                    for name in archive.namelist():
                        if name.startswith(source.name + '/') and name.endswith(('.so', '.dylib')):
                            path = Path(archive.extract(name, root / arch))
                            mac.setdefault(Path(name), []).append((arch, path))
        if not inventory:
            raise ValueError('No wheels compatible with the official Blender runtimes')
        if len({frozenset(abis) for abis in inventory.values()}) > 1:
            raise ValueError(f'Wheel platforms cover different Python ABIs: {inventory}')

        mac_archs = {f.removeprefix('mac-') for f in inventory if f.startswith('mac-')}
        for name, slices in mac.items():
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if {arch for arch, _ in slices} != mac_archs:
                raise ValueError(f'Missing macOS slice for {name}')
            if len(slices) == 1:
                shutil.copyfile(slices[0][1], target)
            else:
                subprocess.run(['lipo', '-create', *(path for _, path in slices), '-output', target],
                               check=True)
            subprocess.run(['lipo', target, '-verify_arch', *(arch for arch, _ in slices)], check=True)

        execute_Export(package, destination, source.name)


def _extract_linux(archive, package):
    # auditwheel places private libraries beside the Python package. Blender installs
    # a single add-on directory, so move those libraries inside it and adjust RPATH.
    for name in archive.namelist():
        path = Path(name)
        if name.endswith('/'):
            continue
        if path.parts[0] == package.name + '.libs':
            target = package / '.libs' / path.name
        elif path.parent == Path(package.name) and path.suffix == '.so':
            target = package / path.name
        else:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(archive.read(name))
        old = subprocess.check_output(['patchelf', '--print-rpath', target], text=True).strip()
        new = old.replace('$ORIGIN/../' + package.name + '.libs',
                          '$ORIGIN' if target.parent != package else '$ORIGIN/.libs')
        if new != old:
            subprocess.run(['patchelf', '--set-rpath', new, target], check=True)

"""Install an official Blender build in the CI job's temporary working directory."""
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import sys
from zipfile import ZipFile


# Checksums from each release's blender-VERSION.sha256.
ARCHIVES = {
    "3.6.7": {
        "linux-x64.tar.xz": "c57587bddfc907e9cc55fba5daf51b6053d117ff0652dc9b0a4504711458238e",
        "macos-arm64.dmg": "8275bffaca91d8c7094f992a47f0e515ce5e19c0a81b66a2914fbfb34e5eebe4",
        "windows-x64.zip": "043053f9c25e9df7ddb94dd2afb06fb6f6cffbb770966d7045247d868cbca46b",
    },
    "4.3.2": {
        "linux-x64.tar.xz": "4da1c956673c0485e63054e563ee69198cc8f80d8157dd7592dffc8a6a5592e6",
        "macos-arm64.dmg": "32297d20fae74683052267b9171945712c0aeb8a499d8e677e44740d3e30a034",
        "windows-x64.zip": "933a62a770c941e316535f893629ff106ddd2aa5a956f0081ec557439073d4f6",
    },
}


def install(version):
    suffix = {("Linux", "x86_64"): "linux-x64.tar.xz", ("Darwin", "arm64"): "macos-arm64.dmg",
              ("Windows", "amd64"): "windows-x64.zip"}[platform.system(), platform.machine().lower()]
    checksum = ARCHIVES[version][suffix]
    archive = Path(f"blender-{version}-{suffix}")
    release = version.rsplit(".", 1)[0]
    url = f"https://download.blender.org/release/Blender{release}/{archive.name}"
    subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--output", archive, url], check=True)
    with archive.open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        if digest.hexdigest() != checksum:
            raise ValueError(f"SHA-256 mismatch for {archive}")
    if sys.platform == "darwin":
        mount = Path.cwd() / "blender-dmg"
        subprocess.run(["hdiutil", "attach", "-readonly", "-nobrowse", "-mountpoint", mount, archive], check=True)
        binary = mount / "Blender.app/Contents/MacOS/Blender"
    elif sys.platform == "win32":
        with ZipFile(archive) as distribution:
            distribution.extractall()
        binary = Path(archive.stem) / "blender.exe"
    else:
        subprocess.run(["tar", "-xf", archive], check=True)
        binary = Path(archive.name.removesuffix(".tar.xz")) / "blender"
    with open(os.environ["GITHUB_ENV"], "a", encoding="utf-8") as env:
        print(f"BLENDER_BIN={binary.resolve(strict=True)}", file=env)


if __name__ == "__main__":
    install(sys.argv[1])

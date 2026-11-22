"""Stage the static site using the sample and links from a published release."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = 'https://github.com/smparker/molecular-blender'


def prepare(tag, destination):
    if not re.fullmatch(r'v[0-9][A-Za-z0-9._-]*', tag):
        raise ValueError(f'Invalid release tag {tag!r}')
    ref = f'refs/tags/{tag}'
    sample = 'water-sto-3g.molden'
    data = subprocess.check_output(['git', 'show', f'{ref}:examples/{sample}'], cwd=ROOT)
    html = (ROOT / 'site/index.html').read_text(encoding='utf-8')
    prefix = f'{REPO}/blob/master/examples/'
    for example in re.findall(re.escape(prefix) + r'([^"]+)', html):
        subprocess.run(['git', 'cat-file', '-e', f'{ref}:examples/{example}'], cwd=ROOT, check=True)
    html = html.replace(prefix, f'{REPO}/blob/{tag}/examples/')
    html = html.replace(f'{REPO}/releases/latest/download/', f'{REPO}/releases/download/{tag}/')
    html = html.replace('Download for Blender</a>', f'Download {tag[1:]} for Blender</a>')

    shutil.copytree(ROOT / 'site', destination)
    (destination / 'index.html').write_text(html, encoding='utf-8')
    (destination / 'samples' / sample).write_bytes(data)
    # Keep external links to the four pages the former Jekyll theme generated.
    for anchor in ('intro', 'why', 'usage', 'gallery'):
        target = f'./#{anchor}'
        (destination / f'{anchor}.html').write_text(
            f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Molecular Blender</title>\n'
            f'<meta http-equiv="refresh" content="0;url={target}"><a href="{target}">Molecular Blender</a>\n',
            encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tag')
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    prepare(args.tag, args.destination)

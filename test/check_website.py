"""Check the built project site, including its images and local navigation."""
import argparse
from html.parser import HTMLParser
from io import BytesIO
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit
from urllib.request import urlopen
from zipfile import ZipFile


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(path.read_text(encoding='utf-8'))
        self.close()

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if 'id' in values:
            anchor = values['id']
            if anchor in self.ids:
                raise ValueError(f'Duplicate HTML id {anchor!r}')
            self.ids.add(anchor)
        if tag == 'img' and 'alt' not in values:
            raise ValueError(f'Missing alt attribute on {values.get("src")}')
        if 'srcset' in values:
            self.links.extend(candidate.split()[0] for candidate in values['srcset'].split(','))
        for key in ('href', 'src'):
            if key in values:
                self.links.append(values[key])


def check_downloads(root, links):
    urls = [link for link in links if link.endswith('/molecular_blender.zip')]
    if len(urls) != 2 or len(set(urls)) != 1:
        raise ValueError('Both add-on download links must point to the same published ZIP')
    match = re.fullmatch(r'https://github\.com/smparker/molecular-blender/releases/download/'
                         r'(v[0-9][A-Za-z0-9._-]*)/molecular_blender\.zip', urls[0])
    if not match:
        raise ValueError(f'Expected a versioned release ZIP, got {urls[0]}')
    tag = match[1]
    examples = [link for link in links if '/examples/' in link]
    if not examples or any(f'/blob/{tag}/examples/' not in link for link in examples):
        raise ValueError(f'Expected example links for {tag}')
    sample = 'water-sto-3g.molden'
    expected = subprocess.check_output(['git', 'show', f'refs/tags/{tag}:examples/{sample}'])
    if (root / 'samples' / sample).read_bytes() != expected:
        raise ValueError(f'Water sample differs from {tag}')
    with urlopen(urls[0], timeout=90) as response, ZipFile(BytesIO(response.read())) as archive:
        if 'molecular_blender/__init__.py' not in archive.namelist() or archive.testzip():
            raise ValueError('The release download must be an intact Blender add-on ZIP')
    print(f'Checked {tag} add-on ZIP, Molden file, and example links')


def check(root, downloads=False):
    root = root.resolve()
    if not (root / 'index.html').is_file():
        raise ValueError('Missing site entry point')
    pages = {path: Page(path) for path in root.rglob('*.html')}
    for path, page in pages.items():
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            if url.path.startswith('/'):
                raise ValueError(f'{path.name}: {link} must work below the project URL')
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            if target.is_dir():
                target /= 'index.html'
            if not target.is_relative_to(root) or not target.is_file():
                raise ValueError(f'{path.name}: missing local target {link}')
            if url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                raise ValueError(f'{path.name}: missing anchor {link}')
    print(f'Checked local links and images in {len(pages)} HTML pages')
    if downloads:
        check_downloads(root, pages[root / 'index.html'].links)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--downloads', action='store_true', help='Verify the staged release against GitHub')
    args = parser.parse_args()
    check(args.root, args.downloads)

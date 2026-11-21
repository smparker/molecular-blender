"""Check the built project site, including its images and local navigation."""
from html.parser import HTMLParser
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit


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
        for key in ('href', 'src'):
            if key in values:
                self.links.append(values[key])


def check(root):
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


if __name__ == '__main__':
    check(Path(sys.argv[1]))

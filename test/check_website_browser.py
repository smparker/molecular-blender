"""Exercise the public site's navigation, layout, and sample download in Chromium."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread

from playwright.sync_api import expect, sync_playwright


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def check_viewport(browser, url, width, height, scale, screenshots):
    page = browser.new_page(viewport={'width': width, 'height': height}, has_touch=width <= 390)
    errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.on('requestfailed', lambda request: errors.append(request.url))
    page.on('response', lambda response: errors.append(response.url) if response.status >= 400 else None)
    page.goto(url)
    if scale != 100:
        page.add_style_tag(content=f':root {{ font-size: {scale}% }}')
    page.evaluate('document.fonts.ready')
    assert page.evaluate('document.fonts.check(\'600 1em "Raleway"\')')
    nav = page.get_by_role('navigation')
    page.keyboard.press('Tab')
    expect(page.get_by_role('link', name='Skip to content')).to_be_focused()
    page.keyboard.press('Enter')
    expect(page.get_by_role('main')).to_be_focused()
    for anchor in ('intro', 'why', 'usage', 'gallery'):
        link = nav.locator(f'[href="#{anchor}"]')
        link.focus()
        page.keyboard.press('Enter')
        expect(page).to_have_url(url + '#' + anchor)
        expect(link).to_have_attribute('aria-current', 'location')
        section = page.locator('#' + anchor)
        heading = section.get_by_role('heading').first
        assert heading.bounding_box()['y'] >= page.evaluate(
            "getComputedStyle(document.querySelector('header')).position === 'sticky' "
            "? document.querySelector('header').getBoundingClientRect().bottom : 0") - 1
        for img in section.locator('img').all():
            img.scroll_into_view_if_needed()
            img.evaluate('(img) => img.decode()')
            bounds = img.bounding_box()
            assert bounds['x'] >= 0 and bounds['x'] + bounds['width'] <= width + 1
        if anchor == 'why' and width in (390, 1280):
            video = section.locator('video')
            video.scroll_into_view_if_needed()
            page.wait_for_function("document.querySelector('video').currentTime > 0.1")
            assert video.evaluate('v => !v.controls && v.loop && v.muted && v.playsInline && v.videoWidth > 0')
            video.evaluate('v => { v.currentTime = v.duration - 0.2; }')
            page.wait_for_function("document.querySelector('video').currentTime < 1")
            toggle = section.get_by_role('button', name='Pause trajectory')
            assert toggle.bounding_box() == video.bounding_box()
            if width == 390:
                toggle.tap()
            else:
                toggle.click()
            expect(video).to_have_js_property('paused', True)
            section.get_by_role('button', name='Play trajectory').click()
            expect(video).to_have_js_property('paused', False)
        assert page.evaluate('document.documentElement.scrollWidth') == width, (width, scale, anchor)
    hero = page.locator('.intro-copy').bounding_box()
    figure = page.locator('.orbital').bounding_box()
    if width == 1280:
        assert figure['x'] >= hero['x'] + hero['width']
    elif width <= 390:
        assert figure['y'] >= hero['y'] + hero['height']
        assert all(box.bounding_box()['height'] >= 44 for box in nav.get_by_role('link').all())
    page.evaluate('scrollTo(0, 0)')
    expect(nav.get_by_role('link', name='Molecular Blender home')).to_have_attribute('aria-current', 'location')
    if screenshots:
        page.screenshot(path=screenshots / f'site-{width}-{scale}.png', full_page=True)
    assert not errors, errors
    page.close()
    print(f'{width} × {height} at {scale}% text passed')


def check_site(root, channel, screenshots):
    with TemporaryDirectory() as tmp:
        # Mount below a project path just like GitHub Pages.
        (Path(tmp) / 'molecular-blender').symlink_to(root.resolve(), target_is_directory=True)
        with ThreadingHTTPServer(('127.0.0.1', 0), partial(QuietHandler, directory=tmp)) as server:
            Thread(target=server.serve_forever, daemon=True).start()
            url = f'http://127.0.0.1:{server.server_port}/molecular-blender/'
            try:
                with sync_playwright() as playwright:
                    browser = playwright.chromium.launch(channel=channel)
                    for dimensions in [(1280, 800, 100), (768, 1024, 100), (390, 844, 100),
                                       (320, 640, 100), (320, 640, 200)]:
                        check_viewport(browser, url, *dimensions, screenshots)
                    page = browser.new_page(reduced_motion='reduce')
                    page.goto(url + '#why')
                    video = page.locator('video')
                    video.scroll_into_view_if_needed()
                    assert video.evaluate('v => v.paused && v.currentTime === 0 && !!v.poster')
                    page.get_by_role('button', name='Play trajectory').focus()
                    page.keyboard.press('Space')
                    page.wait_for_function("document.querySelector('video').currentTime > 0.1")
                    page.keyboard.press('Space')
                    expect(video).to_have_js_property('paused', True)
                    page.close()
                    page = browser.new_page(java_script_enabled=False)
                    page.goto(url)
                    expect(page.get_by_role('link', name='Play the trajectory video')).to_have_attribute('href', 'img/trajectory.mp4')
                    page.get_by_role('link', name='Installation guide').click()
                    expect(page).to_have_url(url + '#usage')
                    with page.expect_download() as result:
                        page.get_by_role('link', name='Download the water sample').click()
                    download = result.value
                    assert download.suggested_filename == 'water-sto-3g.molden'
                    assert download.path().read_bytes() == (root / 'samples' / download.suggested_filename).read_bytes()
                    assert download.path().read_bytes().startswith(b'[Molden Format]')
                    browser.close()
                    print('Installation guide and one-click Molden download passed with JavaScript disabled')
            finally:
                server.shutdown()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--channel', help='Use an installed browser, e.g. chrome')
    parser.add_argument('--screenshots', type=Path, help='Save full-page screenshots into this directory')
    args = parser.parse_args()
    check_site(args.root, args.channel, args.screenshots)

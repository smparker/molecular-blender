// Anchor links work on their own. Highlight the section under the fixed menu.
const links = Array.from(document.querySelectorAll('nav a'));
const sections = links.map(link => document.getElementById(link.hash.slice(1)));

function markCurrentSection() {
  let current = 0;
  sections.forEach((section, index) => {
    const offset = parseFloat(getComputedStyle(section).scrollMarginTop);
    if (section.getBoundingClientRect().top <= offset + 1) current = index;
  });
  links.forEach((link, index) => {
    link.parentElement.classList.toggle('active', index === current);
    if (index === current) link.setAttribute('aria-current', 'location');
    else link.removeAttribute('aria-current');
  });
}

window.addEventListener('scroll', markCurrentSection, { passive: true });
window.addEventListener('resize', markCurrentSection);
markCurrentSection();

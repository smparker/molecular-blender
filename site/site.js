// Cache section positions when the layout changes. Scrolling only looks up an offset.
const header = document.querySelector('header');
const links = [...header.querySelectorAll('a[href^="#"]')];
const sections = links.map(link => document.querySelector(link.hash));
let starts = [];
let lastScroll = Infinity;
let current;

function markSection() {
  const y = window.scrollY;
  // A short last section may never reach the sticky header.
  const index = y > 0 && y >= lastScroll - 1
    ? links.length - 1
    : Math.max(0, starts.findLastIndex(start => start <= y + 1));
  const next = links[index];
  if (current === next) return;
  current?.removeAttribute('aria-current');
  next.setAttribute('aria-current', 'location');
  current = next;
}

function measureLayout() {
  const sticky = getComputedStyle(header).position === 'sticky';
  document.documentElement.style.setProperty('--header-space', (sticky ? header.offsetHeight : 0) + 'px');
  const clearance = parseFloat(getComputedStyle(document.documentElement).scrollPaddingTop);
  starts = sections.map(section => section.getBoundingClientRect().top + window.scrollY - clearance);
  lastScroll = document.documentElement.scrollHeight - window.innerHeight;
  markSection();
}

const layout = new ResizeObserver(measureLayout);
for (const element of [header, ...sections, document.querySelector('footer')]) layout.observe(element);
window.addEventListener('resize', measureLayout);
window.addEventListener('scroll', markSection, { passive: true });
measureLayout();

const trajectory = document.querySelector('#trajectory');
const trajectoryToggle = document.querySelector('#trajectory-toggle');
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');

function playTrajectory() {
  // Browsers can require a tap before starting even muted video.
  trajectory.play().catch(() => {});
}

trajectoryToggle.hidden = false;
trajectoryToggle.addEventListener('click', () => {
  if (trajectory.paused) playTrajectory();
  else trajectory.pause();
});
trajectory.addEventListener('play', () => trajectoryToggle.setAttribute('aria-label', 'Pause trajectory'));
trajectory.addEventListener('pause', () => trajectoryToggle.setAttribute('aria-label', 'Play trajectory'));

function applyMotionPreference() {
  if (reducedMotion.matches) trajectory.pause();
  else playTrajectory();
}

reducedMotion.addEventListener('change', applyMotionPreference);
applyMotionPreference();

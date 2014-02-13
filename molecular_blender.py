import periodictable as pt

import sys

class Atom():
    el = pt.element()
    position = (0.0, 0.0, 0.0)

    def __init__(self, symbol, position):
        self.el = pt.elements[symbol]
        self.position = position

class element():
    vdw = 2.0
    color = (1.0, 1.0, 1.0)
    mass = 1.0
    name = ""
    symbol = ""

    def __init__(self, vdw, col, m, n, s):
        self.vdw = vdw
        self.color = col
        self.mass = m
        self.name= n
        self.symbol = s

elements = {
  'h' : element(1.20, (1.0, 1.0, 1.0), 1.0, "Hydrogen", "h"),
  'c' : element(1.70, (0.2, 0.2, 0.2), 12.0, "Carbon", "c"),
  'n' : element(1.55, (0.0, 0.0, 1.0), 14.0, "Nitrogen", "n"),
  'o' : element(1.52, (1.0, 0.0, 0.0), 16.0, "Oxygen", "o")
}

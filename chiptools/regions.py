from collections import namedtuple
Region = namedtuple("Region", ["start", "end", "direction"])

class Regions:
    def __init__(self, starts, ends, directions=1):
        self.starts = starts
        self.ends = ends
        self.directions = directions

    def __iter__(self):
        return (Region(*t) for t in zip(self.starts, self.ends, self.directions))

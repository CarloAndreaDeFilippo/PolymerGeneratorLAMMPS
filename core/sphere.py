import math as m

class Sphere:

    def __init__(
            self, 
            sigma: float = 1., 
            cm : list[float, float, float] = [0., 0., 0.],
            atomID: int | None = None,
            molID: int | None = None,
            atomType: int | None = None,
            ):
        
        self.sigma = sigma
        self.cm = cm
        self.atomID = atomID
        self.molID = molID
        self.atomType = atomType

    def distance(self, other, Lbox) -> float:

        dist = [0., 0., 0.]

        for ax in range(3):
            dist[ax] = self.cm[ax] - other.cm[ax]
            dist[ax] = dist[ax] - Lbox[ax] * round(dist[ax] / Lbox[ax])
        
        return m.sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2])

    def overlap(self, other, Lbox) -> bool:

        return self.distance(other, Lbox) <= 0.5 * (self.sigma + other.sigma)
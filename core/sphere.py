from enum import Enum, auto
import math as m

class SphereType(Enum):
    GENERIC = auto()
    BEAD = auto()
    PATCH = auto()
    COLLOID = auto()
    SOLVENT = auto()


class Sphere:

    def __init__(self, sigma=1, cm=(0., 0., 0.), sphereType=SphereType.GENERIC):
        
        self.sigma = sigma
        self.cm = cm
        self.sphereType = sphereType

    def distance(self, other, Lbox):

        dist = [0., 0., 0.]

        for ax in range(3):
            dist[ax] = self.cm[ax] - other.cm[ax]
            dist[ax] = dist[ax] - Lbox[ax] * round(dist[ax] / Lbox[ax])
        
        return m.sqrt(dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2])

    def overlap(self, other, Lbox):

        if self.distance(other, Lbox) <= 0.5 * (self.sigma + other.sigma):

            return True

        return False

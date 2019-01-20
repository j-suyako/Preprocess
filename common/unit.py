from lib.simple3D import *
from util.util import transform


class Unit(Polyhedron):

    hole = "void"
    brick = "brick"
    mortar = "mortar"

    def __init__(self, verts, ridge, material):
        super(Unit, self).__init__(verts, ridge)
        self.material = material

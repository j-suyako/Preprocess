from lib.simple3D import *


class Unit(Polyhedron):

    hole = "void"

    def __init__(self, poly=None, points=None, ridge=None, material=None):
        if poly is None:
            if points is not None:
                super(Unit, self).__init__(points, ridge)
            else:
                raise TypeError()
        elif isinstance(poly, Polyhedron):
            self.points = poly.points
            self.ridge = poly.ridge
            self.range = poly.range
        else:
            raise TypeError("")
        self.material = material

    # @property
    # def material(self):
    #     return self.material
    #
    # @material.setter
    # def material(self, material):
    #     self.material = material

    def cut(self, that):
        if isinstance(that, Unit):
            in_unit, out_unit = super(Unit, self).cut(that)
            if in_unit is not None:
                in_unit = Unit(in_unit, self.material)
            if out_unit is not None:
                out_unit = Unit(out_unit, that.material)
            return in_unit, out_unit


if __name__ == "__main__":
    points = np.random.random((10, 3))
    unit = Unit(points)
    assert isinstance(unit, Polyhedron)

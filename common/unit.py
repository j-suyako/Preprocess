from lib.simple3D import *


class Unit(Polyhedron):

    hole = "void"
    brick = "brick"
    mortar = "mortar"

    def __init__(self, index, poly=None, points=None, ridge=None, material=None):
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
        self.index = index
        self.material = material

    # @property
    # def material(self):
    #     return self.material
    #
    # @material.setter
    # def material(self, material):
    #     self.material = material

    def cutby(self, that):
        if isinstance(that, Unit):
            in_unit, out_unit = super(Unit, self).cutby(that)
            if in_unit is not None:
                in_unit = Unit(index=self.index, poly=in_unit, material=self.material)
            if out_unit is not None:
                out_unit = Unit(index=self.index, poly=out_unit, material=that.material)
            return in_unit, out_unit
        else:
            above_poly, down_poly = super(Unit, self).cutby(that)
            above_unit = Unit(index=self.index, poly=above_poly, material=self.material) \
                if above_poly is not None else None
            down_unit = Unit(index=self.index, poly=down_poly, material=self.material)\
                if down_poly is not None else None
            return above_unit, down_unit


if __name__ == "__main__":
    points = np.random.random((10, 3))
    unit = Unit(points)
    assert isinstance(unit, Polyhedron)

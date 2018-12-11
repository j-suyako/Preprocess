from lib.simple3D import *
from util.util import transform


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
            try:
                above_poly, down_poly = super(Unit, self).cutby(that)
                above_unit = Unit(index=self.index, poly=above_poly, material=self.material) \
                    if above_poly is not None else None
                down_unit = Unit(index=self.index, poly=down_poly, material=self.material) \
                    if down_poly is not None else None
                return above_unit, down_unit
            except Exception:
                raise ValueError("points are:\n{}\n\nridge are:\n{}\n\nplane is:\n{}".format(transform(self.points), self.ridge, that))


if __name__ == "__main__":
    points = np.array([[102.81121, 56.7657, 87.45533],
                       [94.25503, 53.92903, 86.51536],
                       [102.00645, 58.62968, 81.12211],
                       [109.83072, 51.60695, 78.89316],
                       [109.45106, 50.38, 76.50364],
                       [106.8728, 56.25869, 80.4479],
                       [105.87575, 56.94286, 79.71741],
                       [106.6709, 52.00454, 74.06481],
                       [107.56315, 50.94866, 74.0888],
                       [96.18606, 51.20256, 88.8235],
                       [103.33226, 51.34123, 90.99724],
                       [105.18233, 48.37969, 90.06792],
                       [105.00402, 48.44222, 90.1998],
                       [104.44012, 49.49415, 90.50178],
                       [100.0126, 51.29744, 72.26513],
                       [96.48603, 51.68801, 72.20569],
                       [92.67727, 52.42107, 81.81529],
                       [92.38668, 50.84859, 79.97423],
                       [95.07553, 53.27744, 73.85617],
                       [95.06625, 53.22555, 73.79419],
                       [101.72806, 57.51207, 86.6481],
                       [101.1781, 58.36262, 84.17336],
                       [99.86463, 57.65028, 85.97272],
                       [104.78916, 48.13627, 90.0],
                       [105.24067, 48.3329, 90.0],
                       [109.22309, 50.0, 76.4786],
                       [106.68935, 50.0, 73.69108],
                       [103.47825, 50.0, 72.88751],
                       [108.17828, 50.0, 83.24655],
                       [96.37461, 50.0, 88.34348],
                       [96.47181, 50.0, 73.24841],
                       [92.6722, 50.0, 79.97992],
                       [101.97973, 50.0, 90.0],
                       [104.48532, 50.0, 90.0]])
    ridge = [[0, 20, 22, 1, 9, 10], [2, 6, 7, 14, 15, 19, 18], [11, 12, 13], [20, 21, 22], [3, 5, 6, 7, 8, 4], [26, 8, 4, 25], [26, 8, 7, 14, 27], [25, 4, 3, 28], [23, 12, 11, 24], [29, 9, 10, 13, 12, 23, 32], [30, 15, 14, 27], [16, 17, 19, 18], [30, 15, 19, 17, 31], [31, 17, 16, 1, 9, 29], [24, 11, 13, 10, 0, 5, 3, 28, 33], [1, 22, 21, 2, 18, 16], [0, 5, 6, 2, 21, 20], [32, 23, 24, 33], [30, 31, 29, 32, 33, 28, 25, 26, 27]]
    unit = Unit(index=0, points=points, ridge=ridge)
    plane = InfinitePlane(normal_vector=np.array([0, 0, 1]), intercept=-90)
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    unit.plot(ax1, alpha=0.2)
    plane.plot(ax1, alpha=0.2)
    above, down = unit.cutby(plane)
    above.plot(ax2)
    down.plot(ax2)
    plt.show()

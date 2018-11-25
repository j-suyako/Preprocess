from lib.simple3D import *


class Cuboid(Polyhedron):
    def __init__(self, origin, side):
        x, y, z = origin
        a, b, c = side
        points = np.array([[px, py, pz] for px in [x, x + a] for py in [y, y + b] for pz in [z, z + c]])
        ridge = [[0, 1, 3, 2], [0, 1, 5, 4], [0, 2, 6, 4], [1, 3, 7, 5], [2, 3, 7, 6], [4, 5, 7, 6]]
        super(Cuboid, self).__init__(points, ridge)
        self.origin = origin
        self.end = (x + a, y + b, z + c)

    def volume(self):
        return np.prod([self.range[i][1] - self.range[i][0] for i in range(3)])

    def centroid(self):
        return (self.range[:, 0] + self.range[:, 1]) / 2

    def contains(self, point):
        around_point = np.around(point, decimals=5)
        return np.all([self.range[i][0] <= around_point[i] <= self.range[i][1] for i in range(3)])

    def cut(self, poly: Polyhedron):
        in_points_index, out_points_index = list(), list()
        for i, point in enumerate(poly.points):
            if self.contains(point):
                in_points_index.append(i)
            else:
                out_points_index.append(i)

        if len(in_points_index) > 0 and len(out_points_index) > 0:
            in_points = poly.points[in_points_index]
            out_points = poly.points[out_points_index]
            intersect_points = np.array([0, 0, 0], ndmin=2)
            for plane in self.planes():
                intersect_point = poly.intersect(plane)[1]
                if len(intersect_point) != 0:
                    intersect_points = np.r_[intersect_points, intersect_point]

            in_points = np.r_[in_points, intersect_points[1:]]
            out_points = np.r_[out_points, intersect_points[1:]]
            try:
                return Polyhedron(in_points), \
                       Polyhedron(out_points)
            except Exception:
                raise ValueError("in points are:\n{}\n\nout points are:\n{}".format(in_points, out_points))

        elif len(in_points_index) > 0:
            return poly, None

        else:  # len(out_points) > 0
            return None, poly

    def __lt__(self, other):
        if not isinstance(other, Cuboid):
            raise TypeError()
        if self.range[2][0] < other.range[2][0]:
            return True
        elif self.range[2][0] > other.range[2][0]:
            return False
        elif self.range[0][0] < other.range[0][0]:
            return True
        elif self.range[0][0] > other.range[0][0]:
            return False
        elif self.range[1][0] < other.range[1][0]:
            return True
        elif self.range[1][0] > other.range[1][0]:
            return False
        return True

    def __repr__(self):
        return "{} -> {}".format(self.origin, self.end)


if __name__ == "__main__":
    hole = Cuboid(origin=(680, 0, 400), side=(640, 240, 500))
    holes = [hole]
    print(holes)
    # poly = Polyhedron(points=np.array([[690.7771291, 109.43022151, 883.6565288],
    #                                    [708.52754338, 0., 883.98987139],
    #                                    [6.88763107e+02, 1.20515338e+02, 1.00000000e+03],
    #                                    [6.29106855e+02, 1.19681577e+02, 8.78883786e+02],
    #                                    [7.08312292e+02, 0.00000000e+00, 1.00000000e+03],
    #                                    [6.33961481e+02, 1.30016296e+02, 1.00000000e+03],
    #                                    [5.71985042e+02, 7.54021411e+01, 8.73937635e+02],
    #                                    [5.60869704e+02, 0.00000000e+00, 1.00000000e+03],
    #                                    [5.64451000e+02, 7.64611761e+01, 1.00000000e+03],
    #                                    [5.68514168e+02, -7.10542736e-15, 8.72926647e+02]]))
    # ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    # ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    # hole.plot(ax1)
    # poly.plot(ax1)
    # in_hole, out_hole = hole.cut(poly)
    # in_hole.plot(ax2)
    # out_hole.plot(ax2)
    # plt.show()

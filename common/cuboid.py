from lib.simple3D import *


class Cuboid(Polyhedron):
    def __init__(self, points, ridge=None):
        super(Cuboid, self).__init__(points, ridge)

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
                    intersect_points = np.r_[intersect_points, poly.intersect(plane)[1]]

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


if __name__ == "__main__":
    poly1 = Cuboid(np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    poly2 = Cuboid(np.array([[-0.5, -0.5, -0.5], [0.5, -0.5, -0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
                             [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5]]))
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    poly1.plot(ax1)
    poly2.plot(ax1)
    first, second = poly1.cut(poly2)
    first.plot(ax2)
    second.plot(ax2)
    plt.show()

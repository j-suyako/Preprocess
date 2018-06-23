import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull, Voronoi, KDTree


def norm(x):
    return x / np.sqrt(np.dot(x, x))


# class Point(object):
#     def __init__(self, *args):
#         if len(args) == 1:
#             try:
#                 self.x, self.y, self.z = args[0]
#             except ValueError:
#                 raise ValueError("only 3d points are accepted")
#         else:
#             try:
#                 self.x, self.y, self.z = args
#             except ValueError:
#                 raise ValueError("only 3d points are accepted")
#         self.vector = np.array([self.x, self.y, self.z])
#
#     def __add__(self, other: np.ndarray):
#         x, y, z = self.vector + other
#         return Point(x, y, z)
#
#     def __sub__(self, other: np.ndarray):
#         x, y, z = self.vector
#         return Point(x, y, z)
#
#     def __str__(self):
#         loc = list(map(str, [self.x, self.y, self.z]))
#         return '(' + ', '.join(loc) + ")"
#
#     def __repr__(self):
#         return str(self)


class Segment(object):

    def __init__(self, start_point, end_point):
        self.start_point = start_point
        self.end_point = end_point
        self.vector = end_point - start_point

    @staticmethod
    def in_same_plane(this, that, eps=1e-5):
        if not isinstance(this, Segment) and not isinstance(that, Segment):
            raise TypeError()
        v = list()
        for first in [this.start_point, this.end_point]:
            for second in [that.start_point, that.end_point]:
                v.append(Segment(first, second))
        for e in v:
            if np.abs(e.vector < eps).all():
                return True
        v1, v2, v3, v4 = v
        return True if (np.abs(np.cross(norm(v1 % v2), norm(v3 % v4))) <= eps).all() else False

    def intersect(self, that):
        if isinstance(that, Segment):
            if not Segment.in_same_plane(self, that):
                return None
            x0, y0, z0 = self.start_point
            x1, y1, z1 = that.start_point
            a0, b0, c0 = self.vector
            a1, b1, c1 = that.vector
            t0 = ((y0 - y1) * a1 - (x0 - x1) * b1) / (a0 * b1 - b0 * a1)
            return np.array([t0 * a0 + x0, t0 * b0 + y0, t0 * c0 + z0]) if 0 <= t0 <= 1 else None
        elif isinstance(that, Plane):
            nor_ver = that.normal_vector
            t1 = (-that.z_ - np.dot(nor_ver, np.array(self.start_point))) / np.dot(nor_ver, nor_ver)
            t2 = (-that.z_ - np.dot(nor_ver, np.array(self.end_point))) / np.dot(nor_ver, nor_ver)
            return self.intersect(Segment(self.start_point + nor_ver * t1, self.end_point + nor_ver * t2))

    def __mod__(self, other):
        if not isinstance(other, Segment):
            raise TypeError()
        return np.cross(self.vector, other.vector)

    def __str__(self):
        return str(self.start_point) + "->" + str(self.end_point)

    def __repr__(self):
        return str(self)


class Plane(object):

    def __init__(self, points):
        """

        :param points: np.ndarray, points must be in same plane
        """
        self.points = points
        n = points.shape[0]
        if n < 3:
            raise ValueError("at least 3 points are needed")
        line1 = points[1] - points[0]
        line2 = points[2] - points[0]
        temp = 3
        while temp < n and ((np.abs(norm(line1) - norm(line2)) < 1e-5).all() or
                            (np.abs(norm(line1) + norm(line2)) < 1e-5).all()):
            line2 = points[temp] - points[0]
            temp += 1
        if temp == n:
            raise ValueError("all points shouldn't be in one line")
        self.normal_vector = norm(np.cross(line1, line2))
        for i in range(3, n):
            if (np.abs(np.dot(self.normal_vector, norm(points[i] - points[0]))) >= 1e-5).any():
                raise ValueError("all points should be in a same plane")
        self.z_ = -np.dot(self.normal_vector, points[0])


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __le__(self, other):
        return self.start < other.start

    def __lt__(self, other):
        return self.start <= other.start

    def __str__(self):
        return "(%.2f, %.2f)" % (self.start, self.end)

    def __repr__(self):
        return str(self)


class Polyhedron(object):
    def __init__(self, points):
        if points.shape[1] != 3:
            raise ValueError("only 3d points are accepted")
        self.points = points
        self.hull = ConvexHull(points)
        self.segments = [Segment(points[simplex1], points[simplex2])
                         for simplex1, simplex2 in zip(self.hull.simplices[:-1], self.hull.simplices[1:])]
        self.x_range = Range(min(points[:, 0]), max(points[:, 0]))
        self.y_range = Range(min(points[:, 1]), max(points[:, 1]))
        self.z_range = Range(min(points[:, 2]), max(points[:, 2]))

    def intersect(self, plane: [np.ndarray, Plane]):
        if isinstance(plane, np.ndarray):
            plane = Plane(plane)
        point0 = plane.points[0]
        dis = list()
        for point in self.points:
            v1, v2 = point - point0, plane.normal_vector
            cos_theta = v1 * v2 / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
            dis.append(np.dot(v1, v1) * cos_theta)
        ranges = [Range(dis[simplex1], dis[simplex2])
                  for simplex1, simplex2 in zip(self.hull.simplices[:-1], self.hull.simplices[1:])]
        intersect_points = list()
        # for r in ranges:
        #     if r.start * r.end <= 0:

        # TODO(suyako): 求出多面体与平面的交点



if __name__ == '__main__':
    p = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]])
    vor = Voronoi(p)
    pass

import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import numpy as np
from scipy.spatial import ConvexHull, Voronoi, KDTree


def norm(x):
    if len(x.shape) == 1:
        return x / np.sqrt(np.dot(x, x))
    else:
        return np.array([norm(e) for e in x])


def angel(vec1, vec2):
    return np.arccos(np.sum(vec1 * vec2) / np.sqrt(np.sum(vec1 * vec1) + np.sum(vec2 * vec2)))


def rotate_z(points, theta):
    """点绕z轴逆时针旋转（从z轴正方向看）theta角

    :param points:
    :param theta:
    :return:
    """
    rotate_matrix = np.asarray([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return np.dot(points, rotate_matrix)


def rotate_x(points, theta):
    rotate_matrix = np.asarray([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
    return np.dot(points, rotate_matrix)


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

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        x0, y0, z0 = self.start_point
        x1, y1, z1 = self.end_point
        fig.plot([x0, x1], [y0, y1], [z0, z1], 'o-')

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
        lines = (points - points[0])[1:]
        line0 = lines[0]
        for line1 in lines[1:]:
            if (np.abs(norm(line0) - norm(line1)) < 1e-5).all() or (np.abs(norm(line0) + norm(line1)) < 1e-5).all():
                continue
            else:
                break
        else:
            raise ValueError("all points shouldn't be in one line")
        # TODO: one line circumstance
        self.normal_vector = norm(np.cross(line0, lines[1]))
        if self.normal_vector[0] < 0:
            self.normal_vector *= -1
        for i in range(3, n):
            if (np.abs(np.dot(self.normal_vector, norm(points[i] - points[0]))) >= 1e-5).any():
                raise ValueError("all points should be in a same plane")
        self.z_ = -np.dot(self.normal_vector, points[0])
        center = np.sum(points, 0) / points.shape[0]
        relative_points = points - center
        theta1 = angel(self.normal_vector[:2], np.array([0, 1]))
        normal_vector_after_rotate = rotate_z(self.normal_vector, theta1)
        theta2 = angel(normal_vector_after_rotate[1:], np.array([0, 1]))
        relative_points = rotate_z(relative_points, theta1)
        relative_points = rotate_x(relative_points, theta2)
        self.hull = ConvexHull(relative_points[:, :2])
        index = self.hull.vertices
        self.segments = [Segment(points[start], points[end]) for start, end in
                         zip(index, np.append(index[1:], index[0]))]

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        index = self.hull.vertices
        index = np.append(index, index[0])
        verts_x_y = self.points[index][:, :2]
        verts_z = self.points[index][:, 2]
        paths = [Path.LINETO] * (len(verts_z) - 2)
        paths.insert(0, Path.MOVETO)
        paths.append(Path.CLOSEPOLY)
        path = Path(verts_x_y, paths)
        patch = patches.PathPatch(path, alpha=0.2)
        fig.add_patch(patch)
        art3d.pathpatch_2d_to_3d(patch, z=verts_z)
        fig.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], marker='^')

    def __lt__(self, other):
        if not isinstance(other, Plane):
            raise TypeError()
        p = (self.normal_vector * 1e5).astype(int)
        p /= 1e5
        q = (other.normal_vector * 1e5).astype(int)
        q /= 1e5
        z1 = int(self.z_ * 1e5) / 1e5
        z2 = int(self.z_ * 1e5) / 1e5
        if p[0] != q[0]:
            return p[0] < q[0]
        if p[1] != q[1]:
            return p[1] < q[1]
        if p[2] != q[2]:
            return p[2] < q[2]
        return z1 < z2

    def __eq__(self, other):
        if not isinstance(other, Plane):
            raise TypeError()
        flag = (np.abs(self.normal_vector - other.normal_vector) < 1e-5).all() and np.abs(self.z_ - other.z_) < 1e-5
        return flag


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
        prob_planes = list()
        for index in self.hull.simplices:
            prob_planes.append(Plane(points[index]))
        # prob_planes = [Plane(points[index]) for index in self.hull.simplices]
        prob_planes.sort()
        self.planes = list()
        for e1, e2 in prob_planes[:-1], prob_planes[1:]:
            if e1 != e2:
                self.planes.append(e1)
        # self.segments = list()
        # for index1, index2, index3 in self.hull.simplices:
        #     self.segments.append(Segment(points[index1], points[index2]))
        #     self.segments.append(Segment(points[index1], points[index3]))
        #     self.segments.append(Segment(points[index2], points[index3]))
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
        intersect_points = np.array([])
        above_plane_points = np.array([])
        down_plane_points = np.array([])
        for s, r in zip(self.segments, ranges):
            for point, dis in zip([s.start_point, s.end_point], [r.start, r.end]):
                if dis > 0:
                    above_plane_points = np.vstack((above_plane_points, point))
                elif dis < 0:
                    down_plane_points = np.vstack((down_plane_points, point))
            if r.start * r.end <= 0:
                intersect_points = np.vstack((intersect_points, s.intersect(plane)))
        return above_plane_points, intersect_points, down_plane_points

    def cutting(self, plane: [np.ndarray, Plane]):
        if isinstance(plane, np.ndarray):
            plane = Plane(plane)
        above, intersect, down = self.intersect(plane)
        return Polyhedron(np.vstack((above, intersect))), Polyhedron(np.vstack((down, intersect)))

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        for plane in self.planes:
            plane.plot(figure=fig)
        # for segment in self.segments:
        #     segment.plot(figure=fig)
        fig.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], marker='x')
        # for point in self.points:
        #     figure.plot(point)


if __name__ == '__main__':
    # p = np.random.rand(10, 3)
    # vor = Voronoi(p)
    # index = 0
    # for i, e in enumerate(vor.regions):
    #     if e and (np.asarray(e) >= 0).all():
    #         index = i
    #         break
    # points = vor.vertices[vor.regions[index]]
    points = np.random.random((15, 3))
    poly = Polyhedron(points)
    poly.plot()
    plt.show()

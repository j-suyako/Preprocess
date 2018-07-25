import matplotlib.patches as patches
from matplotlib.path import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import numpy as np
from scipy.spatial import ConvexHull, Voronoi, Delaunay


INFINITY = np.array([np.inf, np.inf, np.inf])


def norm(x):
    if len(x.shape) == 1:
        return x / np.sqrt(np.dot(x, x))
    else:
        return np.array([norm(e) for e in x])


def angel(vec1, vec2):
    if (np.abs(vec1) < 1e-5).all() or (np.abs(vec1) < 1e-5).all():
        return 0
    return np.arccos(np.sum(vec1 * vec2) / np.sqrt(np.sum(vec1 * vec1) * np.sum(vec2 * vec2)))


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

    def takepoint(self, p, axis: str):
        """取改线段所在直线上的一点"""
        axis = axis.lower()
        axs = ['x', 'y', 'z']
        if axis not in axs:
            raise ValueError()
        index = axs.index(axis)
        return (p - self.start_point[index]) / self.vector[index] * self.vector + self.start_point

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
            if self == that:
                raise ValueError("these two segments are overlapped")
            if not Segment.in_same_plane(self, that):
                return None
            if np.abs(angel(self.vector, that.vector)) < 1e-5:
                return None
            temp = np.r_[self.vector, -that.vector].reshape(2, 3).T
            index = list()
            for e in [[0, 1], [0, 2]]:
                if np.abs(np.linalg.det(temp[e])) > 1e-5:
                    index = e
                    break
            t0 = np.dot(np.linalg.inv(temp[index]), (that.start_point - self.start_point)[index])[0]
            # t0 = ((y0 - y1) * a1 - (x0 - x1) * b1) / (a0 * b1 - b0 * a1)
            return self.vector * t0 + self.start_point if 0 <= t0 <= 1 else None
            # return np.array([t0 * a0 + x0, t0 * b0 + y0, t0 * c0 + z0]) if 0 <= t0 <= 1 else None
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

    def __eq__(self, other):
        if (np.abs(self.start_point - other.start_point) >= 1e-5).any():
            return False
        if (np.abs(self.end_point - other.end_point) >= 1e-5).any():
            return False
        return True

    def __mod__(self, other):
        if not isinstance(other, Segment):
            raise TypeError()
        return np.cross(self.vector, other.vector)

    def __str__(self):
        return str(self.start_point) + "->" + str(self.end_point)

    def __repr__(self):
        return str(self)


class Plane(object):

    def __init__(self, points=None, normal_vector=None, z_=None):
        """

        :param points: np.ndarray, points must be in same plane
        """
        pass
        if points is not None:
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
            for i in range(3):
                if self.normal_vector[i] == 0:
                    continue
                elif self.normal_vector[i] < 0:
                    self.normal_vector *= -1
                    break
                else:
                    break
            for i in range(3, n):
                if (np.abs(np.dot(self.normal_vector, norm(points[i] - points[0]))) >= 1e-5).any():
                    raise ValueError("all points should be in a same plane")
            self.z_ = -np.dot(self.normal_vector, points[0])
            a, b, c = list(map(lambda x: float('%.5f' % x), self.normal_vector))
            self.feature = (a, b, c, float('%.5f' % self.z_))
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
        else:
            if normal_vector is None or z_ is None:
                raise ValueError("a plane is decided by its normal vector and its constant term.")
            self.normal_vector = normal_vector
            self.z_ = z_
            a, b, c = list(map(lambda x: float('%.5f' % x), self.normal_vector))
            self.feature = (a, b, c, float('%.5f' % self.z_))

    def get_loc(self, x=None, y=None, z=None):
        if x is None and y is None and z is None:
            raise ValueError()
        a, b, c = self.normal_vector
        if x is not None:
            if y is not None:
                return (-self.z_ - (a * x + b * y)) / c if c else 1
            elif z is not None:
                return (-self.z_ - (a * x + c * z)) / b if b else 1
        else:
            return (-self.z_ - (b * y + c * z)) / a if a else 1

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        if hasattr(self, 'points'):
            index = self.hull.vertices
            index = np.append(index, index[0])
            verts_x_y = self.points[index][:, :2]
            verts_z = self.points[index][:, 2]
        else:
            xmin, xmax = fig.get_xlim()
            ymin, ymax = fig.get_ylim()
            verts_x_y = np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]])
            a, b, c = self.normal_vector
            verts_z = (-self.z_ - np.dot(verts_x_y, np.array([a, b]))) / c
        paths = [Path.LINETO] * (len(verts_z) - 2)
        paths.insert(0, Path.MOVETO)
        paths.append(Path.CLOSEPOLY)
        path = Path(verts_x_y, paths)
        patch = patches.PathPatch(path, alpha=0.2)
        fig.add_patch(patch)
        art3d.pathpatch_2d_to_3d(patch, z=verts_z)
        return patch
        # fig.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], marker='^')

    def __repr__(self):
        a, b, c, z = list(map(lambda x: '%.3f' % x if x < 0 else ('+%.3f' % x if x > 0 else '+0'),
                              np.append(self.normal_vector, self.z_)))
        a = a[1:] if a[0] == '+' else a
        return "{}x{}y{}z{} = 0".format(a, b, c, z)

    def __lt__(self, other):
        if not isinstance(other, Plane):
            raise TypeError()
        p = (self.normal_vector * 1e5).astype(int)
        q = (other.normal_vector * 1e5).astype(int)
        z1 = int(self.z_ * 1e5)
        z2 = int(self.z_ * 1e5)
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


class FinitePlane(Plane):
    def __init__(self, points):
        super(self)


class SemiInfinitePlane(Plane):
    pass


class InfinitePlane(Plane):
    pass


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
        if points.shape[0] <= 3:
            self.planes = []
        else:
            if points.shape[1] != 3:
                raise ValueError("only 3d points are accepted")
            self.points = points
            self.hull = ConvexHull(points)
            self.volume = self.hull.volume
            tera = Delaunay(self.hull.points)

            plane_feature_dic = dict()  # key值为平面的法向量及常数项构成的集合，value为{ points }中属于该平面上的点的索引
            for index in self.hull.simplices:
                feature = Plane(points[index]).feature
                plane_feature_dic.setdefault(feature, set())
                plane_feature_dic[feature] |= set(index)
            self.planes = list()
            self.points_connection = np.array([[False] * points.shape[0] for _ in range(points.shape[0])])
            for e in plane_feature_dic:
                index = list(plane_feature_dic[e])
                self.planes.append(Plane(points[index]))
                curr_plane = self.planes[-1]
                for e1, e2 in curr_plane.hull.simplices:
                    start, end = (index[e1], index[e2]) if index[e1] < index[e2] else (index[e2], index[e1])
                    self.points_connection[start][end] = True
            self.x_range = Range(min(points[:, 0]), max(points[:, 0]))
            self.y_range = Range(min(points[:, 1]), max(points[:, 1]))
            self.z_range = Range(min(points[:, 2]), max(points[:, 2]))

    def intersect(self, plane: [np.ndarray, Plane]):
        # TODO(suyako): test this func
        if isinstance(plane, np.ndarray):
            plane = Plane(plane)
        x, y = self.points[0][:2]
        point0 = plane.get_loc(x, y)
        # point0 = plane.points[0]

        # 求所有点与平面的相对距离，在法向量方向那一面的点距离为正，另一面为负，若某线段一段距离为正，
        # 另一端距离为负，则该线段与平面必然相交
        dis = list()
        intersect_points = np.array([])
        above_plane_points = np.array([])
        down_plane_points = np.array([])
        for point in self.points:
            v1, v2 = point - point0, plane.normal_vector
            cos_theta = np.sum(v1 * v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
            dis.append(np.dot(v1, v1) * cos_theta)
            if dis[-1] > 0:
                above_plane_points = np.append(above_plane_points, point)
            elif dis[-1] < 0:
                down_plane_points = np.append(down_plane_points, point)
            else:
                intersect_points = np.append(intersect_points, point)
        N = self.points.shape[0]
        for start in range(N):
            for end in range(start + 1, N):
                if self.points_connection[start][end] and dis[start] * dis[end] < 0:
                    intersect_points = np.append(intersect_points,
                                                 Segment(self.points[start], self.points[end]).intersect(plane))
        return above_plane_points.reshape((-1, 3)), \
               intersect_points.reshape((-1, 3)), \
               down_plane_points.reshape((-1, 3))

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
        # fig.scatter(self.points[:, 0], self.points[:, 1], self.points[:, 2], marker='x')
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

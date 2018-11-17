import numpy as np
from scipy.spatial import ConvexHull, Voronoi, Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import matplotlib.patches as patches
from matplotlib.path import Path
from util.math import norm, angel, rotate_x, rotate_z


class Line(object):

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
        if not isinstance(this, Line) and not isinstance(that, Line):
            raise TypeError()
        v = list()
        for first in [this.start_point, this.end_point]:
            for second in [that.start_point, that.end_point]:
                v.append(Line(first, second))
        for e in v:
            if np.abs(e.vector < eps).all():
                return True
        v1, v2, v3, v4 = v
        return True if (np.abs(np.cross(norm(v1 % v2), norm(v3 % v4))) <= eps).all() else False

    def intersect(self, that):
        if isinstance(that, Line):
            if self == that:
                raise ValueError("these two segments are overlapped")
            if not Line.in_same_plane(self, that):
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
            return self.vector * t0 + self.start_point if 0 <= t0 <= 1 else None

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
        if not isinstance(other, Line):
            raise TypeError()
        return np.cross(self.vector, other.vector)

    def __str__(self):
        return str(self.start_point) + "->" + str(self.end_point)

    def __repr__(self):
        return str(self)


class Plane(object):

    def __init__(self, points=None):
        """
        Before construct the plane, this constructor will judge if these points
        1. 3d or not, if not, throw exception;
        2. all in one line or not, if true, throw exception;
        3. all in same plane or not, if not, throw exception;
        :param points: np.ndarray or None
        """
        self.points = points
        if points is None:
            self.normal_vector = None
            self.intercept = None
        else:
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

            normal_vector = norm(np.cross(line0, line1))
            for line in lines:
                if (np.abs(np.dot(normal_vector, norm(line))) >= 1e-4).any():
                    raise ValueError("all points should be in a same plane, input points are:\n{}".format(points))

            normal_vector = np.array([0 if abs(e) < 1e-5 else e for e in normal_vector])
            if list(normal_vector < 0) > list(normal_vector > 0):
                normal_vector *= -1
            self.normal_vector = normal_vector
            self.intercept = -np.dot(self.normal_vector, points[0])

    def takepoint(self, x=None, y=None, z=None):
        if x is None and y is None and z is None:
            raise ValueError()
        a, b, c = self.normal_vector
        if x is not None:
            if y is not None:
                return (-self.intercept - (a * x + b * y)) / c if c else 1
            elif z is not None:
                return (-self.intercept - (a * x + c * z)) / b if b else 1
        else:
            return (-self.intercept - (b * y + c * z)) / a if a else 1

    def parallel_with_axis(self):
        """
        judege if the normal verctor of this plane is parallel with axis
        if true, return the index of corresponding axis, else return None
        :return:
        """
        a, b, c = self.normal_vector
        if a == 1 and b == 0 and c == 0:
            return True, 0
        elif a == 0 and b == 1 and c == 0:
            return True, 1
        elif a == 0 and b == 0 and c == 1:
            return True, 2
        else:
            return False, None

    def dis_to(self, points: np.ndarray):
        """
        get the distance to target points
        :param points: np.ndarray, target points
        :return: distance array
        """
        if len(points.shape) == 1:
            nor_ver = self.normal_vector
            return (-self.intercept - np.dot(nor_ver, points)) / np.dot(nor_ver, nor_ver)
            # x, y = points[:2]
            # point0 = np.array([x, y, self.takepoint(x, y)])
            # v1, v2 = points - point0, self.normal_vector
            # cos_theta = np.sum(v1 * v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
            # return np.sqrt(np.dot(v1, v1)) * cos_theta
        else:
            return np.array([self.dis_to(point) for point in points])

    def contains(self, point, error=1e-5):
        """
        judge if point in this plane
        :param point: np.array, target point
        :param error: float, allowable error
        :return: boolean, true or false
        """
        if np.abs(np.dot(self.normal_vector, point) + self.intercept) > error:
            return False
        return True

    def intersect(self, line: Line):
        """
        get the point where the plane and the line intersect, this problem finally reduced to the problem
        that get the point between two lines
        :param line: Segment, target line
        :return: intersect point
        """
        nor_ver = self.normal_vector
        start_point, end_point = line.start_point, line.end_point
        dis1, dis2 = self.dis_to(start_point), self.dis_to(end_point)
        if abs(dis1) <= 1e-5:
            return start_point
        elif abs(dis2) <= 1e-5:
            return end_point
        elif dis1 * dis2 < 0:
            project_line = Line(start_point + nor_ver * dis1, end_point + nor_ver * dis2)
            return project_line.start_point + project_line.vector * abs(dis1) / (abs(dis1) + abs(dis2))
        else:
            return None

    def plot(self, figure=None):
        pass

    @staticmethod
    def _plot(points, figure: Axes3D):
        x, y, z = points[:, 0], points[:, 1], points[:, 2]
        curr_minx, curr_maxx = figure.get_xlim()
        curr_miny, curr_maxy = figure.get_ylim()
        curr_minz, curr_maxz = figure.get_zlim()
        minx, maxx = min(min(x), curr_minx), max(max(x), curr_maxx)
        miny, maxy = min(min(y), curr_miny), max(max(y), curr_maxy)
        minz, maxz = min(min(z), curr_minz), max(max(z), curr_maxz)
        figure.set_xlim(minx, maxx)
        figure.set_ylim(miny, maxy)
        figure.set_zlim(minz, maxz)
        paths = [Path.LINETO] * (len(z) - 2)
        paths.insert(0, Path.MOVETO)
        paths.append(Path.CLOSEPOLY)
        path = Path(points[:, :2], paths)
        patch = patches.PathPatch(path, alpha=0.2)
        figure.add_patch(patch)
        art3d.pathpatch_2d_to_3d(patch, z=z)

    def __repr__(self):
        nv = np.around(self.normal_vector, 3)
        nv = np.array([e if e != 0 else 0 for e in nv])
        intercept = np.around(self.intercept, 3)
        intercept = intercept if intercept != 0 else 0
        a, b, c, z = list(map(lambda x: '%.3f' % x if x < 0 else ('+%.3f' % x if x > 0 else '+0'),
                              np.append(nv, intercept)))
        a = a[1:] if a[0] == '+' else a
        return "{}x{}y{}z{}=0".format(a, b, c, z)

    def __lt__(self, other):
        if not isinstance(other, Plane):
            raise TypeError()
        p = (self.normal_vector * 1e5).astype(int)
        q = (other.normal_vector * 1e5).astype(int)
        z1 = int(self.intercept * 1e5)
        z2 = int(self.intercept * 1e5)
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
        flag = (np.abs(self.normal_vector - other.normal_vector) < 1e-5).all() and np.abs(self.intercept - other.intercept) < 1e-5
        return flag


class FinitePlane(Plane):
    def __init__(self, points, ridge=None):
        """
        A finite plane is the convexhull defined by the input points, ridge here descript
        the connection relationship between points
        :param points:
        :param ridge:
        """
        super(FinitePlane, self).__init__(points)
        if ridge is None:
            self.ridge = self._get_ridge(points)
        else:
            self.ridge = ridge
            # self.ridge = self._preprocess_ridge(ridge)

    @staticmethod
    def _preprocess_ridge(ridge):
        map = dict()
        for i, r in enumerate(ridge):
            map.setdefault(r[0], i)
        i = 0
        process_ridge = [ridge[i]]
        while ridge[i][1] != ridge[0][0]:
            i = map[ridge[i][1]]
            process_ridge.append(ridge[i])
        return process_ridge

    def _get_ridge(self, points):
        """
        caulate the connection relationship between points
        :param points:
        :return:
        """
        center = np.sum(points, 0) / points.shape[0]
        relative_points = points - center
        theta1 = angel(self.normal_vector[:2], np.array([0, 1]))
        normal_vector_after_rotate = rotate_z(self.normal_vector, theta1)
        theta2 = angel(normal_vector_after_rotate[1:], np.array([0, 1]))
        relative_points = rotate_z(relative_points, theta1)
        relative_points = rotate_x(relative_points, theta2)
        hull = ConvexHull(relative_points[:, :2])
        # res = [[e1, e2] for e1, e2 in zip(hull.vertices[:-1], hull.vertices[1:])]
        # res.append([hull.vertices[-1], hull.vertices[0]])
        return hull.vertices

    def area(self):
        """
        get area of this finite plane
        :return:
        """
        bound_points = self.points[self.ridge]
        edge_vectors = (bound_points - bound_points[0])[1:]

        area = 0
        for i in range(edge_vectors.shape[0] - 1):
            cross = np.cross(edge_vectors[i], edge_vectors[i + 1])
            area += np.abs(np.sqrt(np.dot(cross, cross))) / 2

        return area

    def contains(self, point, error=1e-5):
        """
        judge if point in this plane
        :param point: np.array, target point
        :param error: float, allowable error
        :return: boolean, true or false
        """
        if super(FinitePlane, self).contains(point, error):
            bound_points = self.points[self.ridge]
            edge_vectors = bound_points - point

            area = 0
            for i in range(edge_vectors.shape[0]):
                start, end = i, i + 1 if i < edge_vectors.shape[0] - 1 else 0
                cross = np.cross(edge_vectors[start], edge_vectors[end])
                area += np.abs(np.sqrt(np.dot(cross, cross))) / 2

            if abs(area - self.area()) < 1e-5:
                return True

        return False

    def intersect(self, line: Line):
        intersect_point = super(FinitePlane, self).intersect(line)
        if self.contains(intersect_point):
            return intersect_point
        else:
            return None

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        figure = Axes3D(plt.figure()) if not figure else figure
        index = list(self.ridge)
        index.append(index[0])
        verts = self.points[index]
        super(FinitePlane, self)._plot(verts, figure)


class InfinitePlane(Plane):
    def __init__(self, normal_vector, intercept):
        super(InfinitePlane, self).__init__()
        self.normal_vector = normal_vector
        self.intercept = intercept

    def plot(self, figure=None):
        # TODO(suyako): plot plane that perpendicular to xy plane, not so important
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        figure = Axes3D(plt.figure()) if not figure else figure
        xmin, xmax = figure.get_xlim()
        ymin, ymax = figure.get_ylim()
        verts_x_y = np.array([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax], [xmin, ymin]])
        a, b, c = self.normal_vector
        verts_z = (-self.intercept - np.dot(verts_x_y, np.array([a, b]))) / c
        verts = np.c_[verts_x_y, verts_z]
        return super(InfinitePlane, self)._plot(verts, figure)


class Polyhedron(object):
    def __init__(self, points, ridge=None):
        """

        :param points:
        :param ridge:
        """
        self.points = points
        if ridge is None:
            self.ridge = self._get_ridge(points)
        else:
            self.ridge = ridge
        # self.volume, self.center = Polyhedron.get_centroid(points)
        self.range = np.array([[np.min(points[:, 0]), np.max(points[:, 0])],
                               [np.min(points[:, 1]), np.max(points[:, 1])],
                               [np.min(points[:, 2]), np.max(points[:, 2])]])

    @staticmethod
    def _get_ridge(points):
        if points.shape[1] != 3:
            raise ValueError("only 3d points are accepted")

        hull = None
        try:
            hull = ConvexHull(points)
        except Exception:
            print(points)

        plane_feature_dic = dict()  # key值为平面的法向量及常数项构成的集合，value为{ points }中属于该平面上的点的索引
        for simplice in hull.simplices:
            curr_plane = Plane(points[simplice])
            feature = str(curr_plane)
            plane_feature_dic.setdefault(feature, set())
            plane_feature_dic[feature] |= set(simplice)

        ridge = list()
        for feature in plane_feature_dic:
            simplice = list(plane_feature_dic[feature])
            try:
                curr_plane = FinitePlane(points[simplice])
                ridge.append([simplice[r] for r in curr_plane.ridge])
            except ValueError:
                raise ValueError("unknown error, the plane may be {}, please check the points".format(feature))

        return ridge

        # points_connection = np.array([[0] * points.shape[0] for _ in range(points.shape[0])])
        # for feature in plane_feature_dic:
        #     index = list(plane_feature_dic[feature])
        #     curr_plane = FinitePlane(points[index])
        #     for e1, e2 in curr_plane.ridge:
        #         start, end = (index[e1], index[e2]) if index[e1] < index[e2] else (index[e2], index[e1])
        #         points_connection[start][end] = True
        #
        # return [[i, j] for i in range(points.shape[0]) for j in range(i + 1, points.shape[0])]

    def _decre_ridge(self):
        n = np.size(self.points, 0)
        points_connection = [[False] * n for _ in range(n)]
        for r in self.ridge:
            for i in range(len(r)):
                start, end = i, i + 1 if i < len(r) - 1 else 0
                minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                if not points_connection[minr][maxr]:
                    yield minr, maxr
                    points_connection[minr][maxr] = True

    def volume(self):
        tetras = Delaunay(self.points)  # compute set of elementary tetrahedra
        vt = 0
        for simplic in tetras.simplices:
            tetra = self.points[simplic]
            vol = np.abs(np.linalg.det(tetra[0:3, :] - tetra[3, :]) / 6)
            vt += vol
        return vt

    def centroid(self):
        """
        get the centroid of a 3d convex polyhedron
        :param points: np.ndarray
        :return:
        """
        try:
            tetras = Delaunay(self.points)  # compute set of elementary tetrahedra
            centroid = np.array([0, 0, 0], dtype="float64")
            vt = 0
            for simplic in tetras.simplices:
                tetra = self.points[simplic]
                centi = np.mean(tetra, 0)
                vol = np.abs(np.linalg.det(tetra[0:3, :] - tetra[3, :]) / 6)
                centroid += centi * vol
                vt += vol
            centroid /= vt
            return centroid
        except Exception:
            raise ValueError("check points:\n{}".format(self.points))

    def contains(self, point):
        pass  # TODO

    def planes(self):
        """
        an iterator, for get planes of this polyhedron
        :return:
        """
        for r in self.ridge:
            ridge = list(range(len(r)))
            yield FinitePlane(self.points[r], ridge)

    def _possible_intersect(self, plane: Plane):
        b, corr_axis = plane.parallel_with_axis()
        if b:
            if not self.range[corr_axis][0] < -plane.intercept < self.range[corr_axis][1]:
                return False
            if isinstance(plane, FinitePlane):
                pass  # TODO: return true directly or judge strictly
            else:
                return True
        else:
            pass

    def _intersect_line(self, line: Line):
        pass  # TODO: poly intersect with line

    def _intersect_plane(self, plane: Plane):
        dis_array = plane.dis_to(self.points)
        above_points = self.points[dis_array > 1e-5]
        down_points = self.points[dis_array < -1e-5]

        intersect_points = list()
        for start, end in self._decre_ridge():
            if dis_array[start] * dis_array[end] < 1e-5:
                intersect_point = plane.intersect(Line(self.points[start], self.points[end]))
                if intersect_point is not None:
                    intersect_points.append(list(intersect_point))

        intersect_points = np.array(intersect_points)
        return above_points, intersect_points, down_points

    def _intersect_poly(self, poly):
        pass

    def intersect(self, that):
        if isinstance(that, Line):
            return self._intersect_line(that)
        elif isinstance(that, Plane):
            return self._intersect_plane(that)
        elif isinstance(that, Polyhedron):
            return self._intersect_poly(that)

    def _cutby_plane(self, plane: Plane):
        """
        poly cut by plane
        :param plane:
        :return:
        """
        above_points, intersect_points, down_points = self._intersect_plane(plane)

        above_points = np.vstack((above_points, intersect_points))
        down_points = np.vstack((down_points, intersect_points))

        return Polyhedron(np.array(above_points)), Polyhedron(np.array(down_points))
        # above_map, down_map = dict(), dict()
        # above_count, down_count = 0, 0
        # above_points, down_points = [], []
        # for i, dis in enumerate(dis_array):
        #     if dis >= 0:
        #         above_points.append(self.points[i])
        #         above_map[i] = above_count
        #         above_count += 1
        #     if dis <= 0:
        #         down_points.append(self.points[i])
        #         down_map[i] = down_count
        #         down_count += 1
        # above_ridge, down_ridge = [], []
        # intersect_points = []
        # above_count_flag, down_count_flag = above_count, down_count
        # for r in self.ridge:
        #     start, end = r[0], r[1]
        #     start_dis, end_dis = dis_array[start], dis_array[end]
        #
        #     if start_dis * end_dis < 0:
        #         curr_edge = Segment(self.points[start], self.points[end])
        #         intersect_point = plane.intersect(curr_edge)
        #         intersect_points.append(intersect_point)
        #         above_points.append(intersect_point)
        #         down_points.append(intersect_point)
        #         if start_dis > 0:
        #             above_ridge.append([above_map[start], above_count])
        #             down_ridge.append([down_map[end], down_count])
        #         else:
        #             above_ridge.append([above_map[end], above_count])
        #             down_ridge.append([down_map[start], down_count])
        #         above_count += 1
        #         down_count += 1
        #
        #     elif start_dis < 0:
        #         down_ridge.append([down_map[start], down_map[end]])
        #
        #     elif start_dis > 0:
        #         above_ridge.append([above_map[start], above_map[end]])
        #
        # intersect_plane = FinitePlane(np.array(intersect_points))
        # above_ridge.extend([[e1 + above_count_flag, e2 + above_count_flag] for e1, e2 in intersect_plane.ridge])
        # down_ridge.extend([[e1 + down_count_flag, e2 + down_count_flag] for e1, e2 in intersect_plane.ridge])
        # return Polyhedron(above_points, above_ridge), Polyhedron(down_points, down_ridge)

    def _cutby_poly(self, poly):
        pass  # TODO: poly cut by poly

    def cutby(self, that):
        """
        if this polyhedron could be intersected with this plane, it'll return the two polyhedrons
        else, it'll just return itself
        :param that:
        :return:
        """
        if isinstance(that, Plane):
            return self._cutby_plane(that)
        elif isinstance(that, Polyhedron):
            return self._cutby_poly(that)

    def cut(self, that):
        if isinstance(that, Polyhedron):
            return that._cutby_poly(self)

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        for plane in self.planes():
            plane.plot(figure=fig)


if __name__ == '__main__':
    points = np.random.random((15, 3))
    poly = Polyhedron(points)
    poly.plot()
    plt.show()

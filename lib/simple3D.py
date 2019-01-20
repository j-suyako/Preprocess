from config import LOGGER
import numpy as np
from scipy.spatial import ConvexHull, Voronoi, Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
import matplotlib.patches as patches
from matplotlib.path import Path
from util.math import norm, angel, rotate_x, rotate_z, ridge_renum, rotate
from util.util import transform


class Line(object):

    def __init__(self, start_point: iter, end_point: iter):
        if not isinstance(start_point, np.ndarray):
            start_point = np.array(list(start_point))
        if not isinstance(end_point, np.ndarray):
            end_point = np.array(end_point)
        self.start_point = start_point
        self.end_point = end_point
        self.vector = end_point - start_point

    def contains(self, point):
        vec = point - self.start_point
        t = np.cross(norm(vec), norm(self.vector))
        if np.sqrt(np.dot(t, t)) < 1e-5:
            return True
        return False

    def in_same_plane(self, that, eps=1e-5):
        if not isinstance(that, Line):
            raise TypeError()
        v = list()
        for first in [self.start_point, self.end_point]:
            for second in [that.start_point, that.end_point]:
                v.append(Line(first, second))
        for e in v:
            if np.abs(e.vector < eps).all():
                return True
        v1, v2, v3, v4 = v
        return True if (np.abs(np.cross(norm(v1 % v2), norm(v3 % v4))) <= eps).all() else False
        # v = list()
        # for first in [self.start_point, self.end_point]:
        #     for second in [that.start_point, that.end_point]:
        #         v.append(Segment(first, second))
        # for e in v:
        #     if np.abs(e.vector < eps).all():
        #         return True
        # v1, v2, v3, v4 = v
        # return True if (np.abs(np.cross(norm(v1 % v2), norm(v3 % v4))) <= eps).all() else False

    def intersect(self, that):
        if isinstance(that, Line):
            if self == that:
                raise ValueError("these two segments are overlapped")
            # if not self.in_same_plane(that):
            #     return None
            if (np.abs(np.cross(norm(self.vector), norm(that.vector))) <= 1e-5).all():
                return None
            temp = np.r_[self.vector, -that.vector].reshape(2, 3).T
            index = list()
            for e in [[0, 1], [0, 2], [1, 2]]:
                if np.abs(np.linalg.det(temp[e])) > 1e-5:
                    index = e
                    break
            try:
                t0, t1 = np.dot(np.linalg.inv(temp[index]), (that.start_point - self.start_point)[index])
                flag1 = 0 < t0 < 1 if isinstance(self, Segment) else t0 > 0
                flag2 = 0 < t1 < 1 if isinstance(that, Segment) else t1 > 0
                if flag1 and flag2:
                    return self.vector * t0 + self.start_point
                return None
            except np.linalg.LinAlgError:
                a = 1

    def plot(self, figure=None):
        pass

    def _plot(self, figure, start, end):
        x0, y0, z0 = start
        x1, y1, z1 = end
        figure.plot([x0, x1], [y0, y1], [z0, z1], 'o-')

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
        return type(self).__name__ + ": " + str(self.start_point) + "->" + str(self.end_point)

    def __repr__(self):
        return str(self)


class Segment(Line):

    def __init__(self, start_point: iter, end_point: iter):
        super(Segment, self).__init__(start_point=start_point, end_point=end_point)

    def takepoint(self, p, axis: str):
        """取该线段所在直线上的一点"""
        axis = axis.lower()
        axs = ['x', 'y', 'z']
        if axis not in axs:
            raise ValueError()
        index = axs.index(axis)
        return (p - self.start_point[index]) / self.vector[index] * self.vector + self.start_point

    # def in_same_plane(self, that, eps=1e-5):
    #     if not isinstance(self, Segment) and not isinstance(that, Segment):
    #         raise TypeError()
    #     v = list()
    #     for first in [self.start_point, self.end_point]:
    #         for second in [that.start_point, that.end_point]:
    #             v.append(Segment(first, second))
    #     for e in v:
    #         if np.abs(e.vector < eps).all():
    #             return True
    #     v1, v2, v3, v4 = v
    #     return True if (np.abs(np.cross(norm(v1 % v2), norm(v3 % v4))) <= eps).all() else False

    def contains(self, point):
        """

        :param point:
        :return:
        """
        if super(Segment, self).contains(point):
            vec0 = point - self.start_point
            vec1 = point - self.end_point
            if np.dot(vec1, self.vector) < 0 < np.dot(vec0, self.vector):
                return True
        return False

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        super(Segment, self)._plot(fig, self.start_point, self.end_point)


class HalfLine(Line):

    def __init__(self, start_point: iter, end_point: iter):
        super(HalfLine, self).__init__(start_point=start_point, end_point=end_point)

    def contains(self, point: iter):
        if super(HalfLine, self).contains(point):
            vec0 = point - self.start_point
            if np.dot(vec0, self.vector) > 0:
                return True
        return False

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        super(HalfLine, self)._plot(fig, self.start_point, self.end_point)


class Plane(object):

    def __init__(self, verts=None, decimals=5):
        """
        Before construct the plane, this constructor will judge if these points
        1. 3d or not, if not, throw exception;
        2. all in one line or not, if true, throw exception;
        3. all in same plane or not, if not, throw exception;
        :param verts: np.ndarray or None
        """
        self.verts = verts
        if verts is None:
            self.normal_vector = None
            self.intercept = None
        else:
            n = verts.shape[0]
            if n < 3:
                raise ValueError("at least 3 points are needed")

            self.verts = np.around(verts, decimals=decimals)
            lines = (self.verts - self.verts[0])[1:]
            line0 = lines[0]
            # for line1 in lines[1:]:
            #     if (np.abs(norm(line0) - norm(line1)) < 1e-5).all() or (np.abs(norm(line0) + norm(line1)) < 1e-5).all():
            #         continue
            #     else:
            #         break
            # else:
            #     raise ValueError("all points shouldn't be in one line, input points are:\n{}".format(self.points))

            normal_vector = norm(np.cross(line0, lines[1]))
            # for line in lines:
            #     if (np.around(np.abs(np.dot(normal_vector, norm(line))), 5) > 1e-5).any():
            #         raise ValueError("all points should be in a same plane, input points are:\n{}".format(self.points))

            normal_vector = np.array([0 if abs(e) < 1e-5 else e for e in normal_vector])
            if list(normal_vector < 0) > list(normal_vector > 0):
                normal_vector *= -1
            self.normal_vector = normal_vector
            self.intercept = -np.dot(self.normal_vector, verts[0])
        self.range = np.array([[np.float('-inf'), np.float('inf')],
                               [np.float('-inf'), np.float('inf')],
                               [np.float('-inf'), np.float('inf')]])

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

    def project(self, points: np.ndarray):
        if len(points.shape) == 1:
            dis = self.dis_to(points)
            return self.normal_vector * dis + points
        else:
            return np.array([list(self.project(point)) for point in points])

    def contains(self, points, error=1e-5):
        """judge if point in this plane

        :param points: np.array, target point
        :param error: float, allowable error
        :return: boolean, true or false
        """
        if np.abs(np.dot(self.normal_vector, points) + self.intercept) > error:
            return False
        return True

    def intersect(self, line: Segment):
        """
        get the point where the plane and the line intersect, this problem finally reduced to the problem
        that get the point between two lines
        :param line: Segment, target line
        :return: intersect point
        """
        nor_ver = self.normal_vector
        start_point, end_point = line.start_point, line.end_point
        dis1, dis2 = np.around([self.dis_to(start_point), self.dis_to(end_point)], 5)
        if dis1 == 0 and dis2 == 0:
            return None
        elif dis1 == 0:
            return start_point
        elif dis2 == 0:
            return end_point
        elif dis1 * dis2 < 0:
            project_line = Segment(start_point + nor_ver * dis1, end_point + nor_ver * dis2)
            return project_line.start_point + project_line.vector * abs(dis1) / (abs(dis1) + abs(dis2))
        else:
            return None

    def plot(self, figure=None, **kwargs):
        pass

    @staticmethod
    def _plot(points, figure: Axes3D, **kwargs):
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
        patch = patches.PathPatch(path, **kwargs)
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
    def __init__(self, verts, ridge=None):
        """A finite plane is the convexhull defined by the input points, ridge here descript the connection relationship between points

        :param verts:
        :param ridge:
        """
        super(FinitePlane, self).__init__(verts)
        if ridge is None:
            self.ridge = self._get_ridge(verts)
        else:
            self.ridge = ridge
            # self.ridge = self._preprocess_ridge(ridge)
        self.range = np.array([[np.min(verts[:, 0]), np.max(verts[:, 0])],
                               [np.min(verts[:, 1]), np.max(verts[:, 1])],
                               [np.min(verts[:, 2]), np.max(verts[:, 2])]])

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
        try:
            hull = ConvexHull(relative_points[:, :2])
        except Exception:
            raise ValueError("input points are:\n{}".format(points))
        # res = [[e1, e2] for e1, e2 in zip(hull.vertices[:-1], hull.vertices[1:])]
        # res.append([hull.vertices[-1], hull.vertices[0]])
        return hull.vertices

    def lines(self):
        for i in range(len(self.ridge)):
            start, end = self.ridge[i], self.ridge[(i + 1) if i < len(self.ridge) - 1 else 0]
            yield Segment(start_point=self.verts[start], end_point=self.verts[end])

    def area(self):
        """
        get area of this finite plane
        :return:
        """
        bound_points = self.verts[self.ridge]
        edge_vectors = (bound_points - bound_points[0])[1:]

        area = 0
        for i in range(edge_vectors.shape[0] - 1):
            cross = np.cross(edge_vectors[i], edge_vectors[i + 1]).astype(float)
            area += np.abs(np.sqrt(np.dot(cross, cross))) / 2

        return area

    def contains(self, points, error=1e-5):
        """
        judge if point in this plane
        :param points: np.array, target point
        :param error: float, allowable error
        :return: boolean, true or false
        """
        if len(points.shape) == 1:
            if super(FinitePlane, self).contains(points):
                midpoint = (self.verts[self.ridge[0]] + self.verts[self.ridge[1]]) / 2
                half_line = HalfLine(start_point=points, end_point=midpoint)
                count = 0
                for line in self.lines():
                    if line.contains(points):
                        return True
                    count += 1 if half_line.intersect(line) is not None else 0

                return True if (count & 1) == 1 else False
            return False
        else:
            return [self.contains(point) for point in points]

    def intersect(self, line: Segment):
        intersect_point = super(FinitePlane, self).intersect(line)
        if intersect_point is not None and self.contains(intersect_point):
            return intersect_point
        else:
            return None

    def plot(self, figure=None, **kwargs):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        figure = Axes3D(plt.figure()) if not figure else figure
        index = list(self.ridge)
        index.append(index[0])
        verts = self.verts[index]
        super(FinitePlane, self)._plot(verts, figure, **kwargs)


class InfinitePlane(Plane):
    def __init__(self, normal_vector, intercept):
        super(InfinitePlane, self).__init__()
        self.normal_vector = normal_vector
        self.intercept = intercept

    def plot(self, figure=None, **kwargs):
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
        return super(InfinitePlane, self)._plot(verts, figure, **kwargs)


class Polyhedron(object):
    def __init__(self, verts, ridge=None, decimals=8):
        """

        :param verts:
        :param ridge:
        """
        self.verts = np.around(verts, decimals=decimals)
        if ridge is None:
            self.ridge = self._get_ridge()
        else:
            self.ridge = ridge
        # self.volume, self.center = Polyhedron.get_centroid(points)
        self.range = np.c_[np.min(verts, 0), np.max(verts, 0)]

    def _get_ridge(self):
        if self.verts.shape[1] != 3:
            raise ValueError("only 3d points are accepted")

        hull = None
        try:
            hull = ConvexHull(self.verts)
        except Exception:
            print(self.verts)

        plane_feature_dic = dict()  # key值为平面的法向量及常数项构成的集合，value为{ points }中属于该平面上的点的索引
        for simplice in hull.simplices:
            curr_plane = Plane(self.verts[simplice])
            feature = str(curr_plane)
            plane_feature_dic.setdefault(feature, set())
            plane_feature_dic[feature] |= set(simplice)

        ridge = list()
        for feature in plane_feature_dic:
            simplice = list(plane_feature_dic[feature])
            try:
                curr_plane = FinitePlane(self.verts[simplice])
                ridge.append([simplice[r] for r in curr_plane.ridge])
            except ValueError:
                raise ValueError("unknown error, the plane may be {}, please check the points".format(feature))

        return ridge

    def _decre_ridge(self):
        n = np.size(self.verts, 0)
        points_connection = [[False] * n for _ in range(n)]
        for r in self.ridge:
            for i in range(len(r)):
                start, end = i, i + 1 if i < len(r) - 1 else 0
                minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                if not points_connection[minr][maxr]:
                    yield minr, maxr
                    points_connection[minr][maxr] = True

    def volume(self):
        tetras = Delaunay(self.verts)  # compute set of elementary tetrahedra
        vt = 0
        for simplic in tetras.simplices:
            tetra = self.verts[simplic]
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
            tetras = Delaunay(self.verts)  # compute set of elementary tetrahedra
            centroid = np.array([0, 0, 0], dtype="float64")
            vt = 0
            for simplic in tetras.simplices:
                tetra = self.verts[simplic]
                centi = np.mean(tetra, 0)
                vol = np.abs(np.linalg.det(tetra[0:3, :] - tetra[3, :]) / 6)
                centroid += centi * vol
                vt += vol
            centroid /= vt
            return centroid
        except Exception:
            raise ValueError("check points:\n{}".format(self.verts))

    def contains(self, point):
        pass  # TODO

    def pan(self, vector: iter):
        if not isinstance(vector, np.ndarray):
            vector = np.array(list(vector), dtype=np.float)
        self.verts += vector
        self.range = (self.range.T + vector).T
        return self

    def planes(self):
        """
        an iterator, for get planes of this polyhedron
        :return:
        """
        for r in self.ridge:
            ridge = list(range(len(r)))
            yield FinitePlane(self.verts[r], ridge)

    def rotate(self, centroid, theta):
        self.verts = rotate(center=centroid, points=self.verts, theta=theta)
        self.range = np.c_[np.min(self.verts, 0), np.max(self.verts, 0)]

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

    def _intersect_line(self, line: Segment, return_key=False):
        points = list()
        key = list()
        for r in self.ridge:
            plane = FinitePlane(self.verts[r], list(range(len(r))))
            point = plane.intersect(line=line)
            if point is not None:
                points.append(list(point))
                key.append(tuple(r))
        if return_key:
            return np.array(points), key
        return np.array(points)

    def _intersect_plane(self, plane: Plane, dis_array=None, return_key=False):
        if dis_array is None:
            dis_array = np.around(plane.dis_to(self.verts), 5)
        if np.all(dis_array >= 0) or np.all(dis_array <= 0):
            return np.array([])
        # intersect_points = [list(e) for e in self.points[dis_array == 0] if plane.contains(e)]
        intersect_points = list()
        key = list()
        for start, end in self._decre_ridge():
            if dis_array[start] * dis_array[end] < 0:
                vector = self.verts[end] - self.verts[start]
                intersect_point = self.verts[start] + vector * abs(dis_array[start]) / \
                                  (abs(dis_array[start]) + abs(dis_array[end]))
                if plane.contains(intersect_point):
                    intersect_points.append(list(intersect_point))
                    key.append((start, end))

        intersect_points = np.array(intersect_points)
        if return_key:
            return intersect_points, key
        return intersect_points

    def _intersect_poly(self, poly):
        pass

    def intersect(self, that, return_key=False):
        if isinstance(that, Segment):
            return self._intersect_line(that, return_key=return_key)
        elif isinstance(that, Plane):
            return self._intersect_plane(that, return_key=return_key)
        elif isinstance(that, Polyhedron):
            return self._intersect_poly(that)

    def _cutby_plane(self, plane: Plane):
        """poly cut by plane

        :param plane:
        :return:
        """
        dis_array = np.around(plane.dis_to(self.verts), 5)
        if np.all(dis_array >= 0):
            return self, None
        elif np.all(dis_array <= 0):
            return None, self
        else:
            above_points = self.verts[dis_array > 0]
            down_points = self.verts[dis_array < 0]
            intersect_points, key = self._intersect_plane(plane, dis_array=dis_array, return_key=True)

            n = self.verts.shape[0]
            mapping = dict(zip(key, list(range(n, n + len(key)))))
            # try:
            #     index = 0
            #     for start, end in self._decre_ridge():
            #         if index == m:
            #             break
            #         if Line(self.points[start], self.points[end]).contains(intersect_points[index]):
            #             mapping[(start, end)] = n + index
            #             index += 1
            # except IndexError:
            #     raise ValueError()

            ridge = [[], []]
            for r in self.ridge:
                flag = list()
                new_r = list()
                for i in range(len(r)):
                    start, end = i, i + 1 if i < len(r) - 1 else 0
                    minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                    new_r.append(r[start])
                    if (minr, maxr) in mapping:
                        flag.append(start)
                        new_r.append(mapping[(minr, maxr)])
                if len(flag) < 2:
                    ridge[0 if dis_array[r[0]] > 0 else 1].append(r)
                else:
                    r0, r1 = new_r[:flag[0] + 2] + new_r[flag[1] + 2:], new_r[flag[0] + 1:flag[1] + 3]
                    if dis_array[r[0]] < 0:
                        r0, r1 = r1, r0
                    ridge[0].append(r0)
                    ridge[1].append(r1)

            r = [e + n for e in FinitePlane(intersect_points).ridge]
            ridge[0].append(r)
            ridge[1].append(r)

            above_points = np.vstack((above_points, intersect_points))
            down_points = np.vstack((down_points, intersect_points))

            # r.sort()
            above_index = list(np.arange(n)[dis_array > 0]) + list(range(n, n + len(r)))
            down_index = list(np.arange(n)[dis_array < 0]) + list(range(n, n + len(r)))
            above_poly = Polyhedron(above_points, ridge=ridge_renum(above_index, ridge=ridge[0]))
            down_poly = Polyhedron(down_points, ridge=ridge_renum(down_index, ridge=ridge[1]))
            return above_poly, down_poly

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

    def plot(self, figure=None, **kwargs):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        fig = Axes3D(plt.figure()) if not figure else figure
        for plane in self.planes():
            plane.plot(figure=fig, **kwargs)

    def to_stream(self):
        lines = list()
        lines.append(self.centroid())
        lines.append(self.verts.shape[0])
        lines += list(self.verts)
        lines.append(len(self.ridge))
        lines += list(self.ridge)
        return lines


if __name__ == '__main__':
    points = np.array([[50.93043, 56.38548, 69.67531],
                       [52.89435, 56.47037, 66.84145],
                       [58.01927, 57.85984, 66.73233],
                       [67.65346, 64.10058, 89.16457],
                       [67.38912, 64.16193, 90.0],
                       [54.22034, 60.54667, 90.0]])
    ridge = [0, 1, 2, 3, 4, 5]
    plane = FinitePlane(verts=points, ridge=ridge)
    line = Segment(start_point=(50, 50, 0), end_point=(50, 50, 90))
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    attribute = {"alpha": 0.2}
    # poly.plot(ax1, **attribute)
    plane.plot(ax1, **attribute)
    line.plot(ax1)
    intersect = plane.intersect(line)
    if intersect is not None:
        halfline = HalfLine(start_point=intersect, end_point=(points[0]+points[1])/2)
        halfline.plot(ax1)
    # above, down = poly.cutby(plane)
    # above.plot(ax2, alpha=0.2)
    # down.plot(ax2, alpha=0.2)
    plt.show()

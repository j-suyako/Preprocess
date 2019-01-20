from lib.simple3D import *
import math
from util.math import rotate
from typing import List


OPTIONAL_DIRECTION = ["pos_xy", "neg_xy", "pos_xz", "neg_xz", "pos_yz", "neg_yz"]


def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.intercept) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


class Cuboid(Polyhedron):
    ridge = [[0, 1, 3, 2], [0, 1, 5, 4], [0, 2, 6, 4], [1, 3, 7, 5], [2, 3, 7, 6], [4, 5, 7, 6]]

    def clone(self):
        clone = Cuboid(self.origin, self.side, self.ele_size)
        if hasattr(self, "points"):
            clone.points = self.points
            clone._update()
        return clone

    def __init__(self, origin, side, ele_size=5):
        if not isinstance(origin, np.ndarray):
            origin = np.array(list(origin), dtype=np.float)
        if not isinstance(side, np.ndarray):
            side = np.array(list(side), dtype=np.float)
        x, y, z = origin
        a, b, c = side
        verts = np.array([[px, py, pz] for px in [x, x + a] for py in [y, y + b] for pz in [z, z + c]], dtype=np.float)
        super(Cuboid, self).__init__(verts, Cuboid.ridge)
        self.origin = origin
        self.side = side
        self.ele_size = ele_size

    def _points_initial(self):
        def func(m, index, d):
            if index == 0:
                p = m + (index + 0.5 * np.random.random()) * d
            elif i == 1:
                p = m + (index + 0.5 + 0.5 * np.random.random()) * d
            elif i == nx - 2:
                p = m + (index + 0.5 * np.random.random()) * d
            elif i == nx - 1:
                p = m + (index + 0.5 + 0.5 * np.random.random()) * d
            else:
                p = m + (index + np.random.random()) * d
            return p

        mx, my, mz = self.origin
        nx, ny, nz = [math.ceil(e / self.ele_size) for e in self.side]
        # nx, ny, nz = np.ceil(self.side / self.ele_size)
        dx, dy, dz = [e[0] / e[1] for e in zip(self.side, [nx, ny, nz])]
        bound_points = list()
        for k in range(nz):
            if nz == 3 and k == 1:
                continue
            for j in range(ny):
                if ny == 3 and j == 1:
                    continue
                for i in range(nx):
                    if nx == 3 and i == 1:
                        continue
                    # if i != 0 and i != nx - 1 and j != 0 and j != ny - 1 and k != 0 and k != nz - 1:
                    #     continue
                    # else:
                    points_x = mx + (i + np.random.random()) * dx if nx == 3 else func(mx, i, dx)
                    points_y = my + (j + np.random.random()) * dy if ny == 3 else func(my, j, dy)
                    points_z = mz + (k + np.random.random()) * dz if nz == 3 else func(mz, k, dz)
                    bound_points.append([points_x, points_y, points_z])

        return np.array(bound_points)

    def initial(self, points=None):
        self.points = self._points_initial() if points is None else points
        self._update()
        return self

    def volume(self):
        return np.prod(self.side)

    def centroid(self):
        return self.origin + self.side / 2

    def contains(self, points):
        if len(points.shape) == 1:
            around_point = np.around(points, decimals=5)
            return np.all([self.range[i][0] <= around_point[i] <= self.range[i][1] for i in range(3)])
        else:
            return [self.contains(point) for point in points]

    def pan(self, vector: iter):
        super(Cuboid, self).pan(vector)
        self.origin += vector
        if hasattr(self, "points"):
            self.points += vector
        return self

    def rotate(self, centroid, direction: str=None):
        if direction is None:
            return
        if direction not in OPTIONAL_DIRECTION:
            raise ValueError()

        # get rotate theta
        index = OPTIONAL_DIRECTION.index(direction)
        theta = [0] * 3
        theta[index // 2] = np.pi / 2 * (1 if index & 1 == 0 else -1)

        # rotate verts && points && swap
        super(Cuboid, self).rotate(centroid, theta)
        if hasattr(self, "points"):
            self.points = rotate(centroid, self.points, theta)
            self._swap(direction)

        self.origin = np.min(self.verts, 0)
        self.side = np.max(self.verts, 0) - np.min(self.verts, 0)

        # return self

    def _swap(self, direction):
        if direction not in OPTIONAL_DIRECTION:
            raise ValueError()

        if direction == "pos_xy":
            self.left_index, self.right_index, self.front_index, self.back_index = \
                self.back_index, self.front_index, self.left_index, self.right_index
        elif direction == "neg_xy":
            self.left_index, self.right_index, self.front_index, self.back_index = \
                self.front_index, self.back_index, self.right_index, self.left_index
        elif direction == "pos_xz":
            self.up_index, self.down_index, self.left_index, self.right_index = \
                self.right_index, self.left_index, self.up_index, self.down_index
        elif direction == "neg_xz":
            self.up_index, self.down_index, self.left_index, self.right_index = \
                self.left_index, self.right_index, self.down_index, self.up_index
        elif direction == "pos_yz":
            self.up_index, self.down_index, self.front_index, self.back_index = \
                self.back_index, self.front_index, self.up_index, self.down_index
        elif direction == "neg_yz":
            self.up_index, self.down_index, self.front_index, self.back_index = \
                self.front_index, self.back_index, self.down_index, self.up_index

    def assemble(self, that):
        if not isinstance(that, Cuboid):
            raise TypeError()
        if not hasattr(self, "points"):
            raise TypeError("Must build cuboid {} first.".format(self))
        if not hasattr(that, "points"):
            raise TypeError("Must build cuboid {} first.".format(that))

        # touching = self.is_touching(that)
        # if not touching:
        #     raise ValueError("Check cuboids {}, {}, it seems that they are not aligning or not touching properly.".
        #                      format(self, that))
        all_verts = np.r_[self.verts, that.verts]
        origin = np.min(all_verts, 0)
        side = np.max(all_verts, 0) - origin
        res = Cuboid(origin=origin, side=side)
        res.points = np.r_[self.points, that.points]
        res._update()
        return res

    def is_touching(self, that):
        relative_direction = norm(that.centroid() - self.centroid())
        if np.all(np.abs(relative_direction) - 1):
            return False
        j = np.argwhere(relative_direction).flatten()[0]
        if np.abs(that.centroid() - self.centroid())[j] != (self.side + that.side)[j] / 2:
            return False

        return True

    def clarify_bound_with(self, that):
        relative_direction = norm(that.centroid() - self.centroid())
        j = np.argwhere(relative_direction).flatten()[0]
        index = [["right_index", "back_index", "up_index"],
                 ["left_index", "front_index", "down_index"]]
        intercepts = [-self.origin - self.side, self.origin]
        i = 1 if sum(relative_direction) == -1 else 0
        touching_face = InfinitePlane(normal_vector=relative_direction, intercept=intercepts[i][j])
        substitute = that.points[getattr(that, index[1 - i][j])]
        self.points[getattr(self, index[i][j])] = mirror_projection(substitute, touching_face)
        # self._update()

    def _update(self):
        nx, ny, nz = np.ceil(self.side / self.ele_size)
        dx, dy, dz = [e[0] / e[1] for e in zip(self.side, [nx, ny, nz])]
        temp = (self.points - self.origin) / np.array([dx, dy, dz])
        self.left_index = np.argwhere(temp[:, 0] <= 1).flatten()
        self.right_index = np.argwhere(temp[:, 0] >= nx - 1).flatten()
        self.front_index = np.argwhere(temp[:, 1] <= 1).flatten()
        self.back_index = np.argwhere(temp[:, 1] >= ny - 1).flatten()
        self.down_index = np.argwhere(temp[:, 2] <= 1).flatten()
        self.up_index = np.argwhere(temp[:, 2] >= nz - 1).flatten()
        index = [self.left_index, self.right_index, self.front_index, self.back_index, self.down_index, self.up_index]
        self.bound_index = np.unique(np.concatenate(index))

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
        return "{}:{}->{}".format(type(self).__name__, self.origin, self.origin + self.side)


if __name__ == "__main__":
    hole = Cuboid(origin=(15, 15, 0), side=(35, 35, 90))
    brick = Cuboid(origin=(50, 15, 0), side=(35, 35, 90))
    hole.initial()
    brick.initial()
    brick.clarify_bound_with(hole)
    hole.assemble(brick)

from common import Unit, Cuboid
from lib.simple3D import *
from typing import List
from util.logging import logging


# TODO: modify the initiation method of points in model
def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.intercept) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


class Block(Cuboid):
    brick = "brick"
    morta = "morta"
    attribute = {
        brick: {"alpha": 0.2, "facecolor": (0.42, 0.1, 0.05)},
        morta: {"alpha": 0.5, "facecolor": (0.75, 0.75, 0.75)},
    }

    def __init__(self, bound: Cuboid, holes: List[Cuboid]=None, ele_size: float=5, material: str="brick"):
        """

        :param bound:
        :param holes:
        :param ele_size:
        :param material:
        """
        if material not in [Block.brick, Block.morta]:
            raise ValueError("Only brick and morta are accepted, your input is {} here.".format(material))
        super(Block, self).__init__(origin=bound.origin, side=bound.side, ele_size=ele_size)
        self.holes = holes
        self.material = material

    def initial(self, points=None, separation=True):
        if separation:
            super(Block, self).initial(points)
        else:
            self.points = super(Block, self)._points_initial()
            self._update()
        # LOGGER.info("{} has {} points.".format(self, self.points.shape[0]))
        # return self

    # @logging("dividing")
    def _divide(self):
        # LOGGER.info("start dividing...")
        index = np.argwhere(self.holes[0].side == self.side)[0][0]
        diagonal_point = self.origin + self.side
        p, q = {0, 1, 2} - {index}
        x, y = [self.origin[p], diagonal_point[p]], [self.origin[q], diagonal_point[q]]
        hole_set = set()
        for hole in self.holes:
            minx, miny = hole.origin[[p, q]]
            maxx, maxy = (hole.origin + hole.side)[[p, q]]
            x.extend([minx, maxx])
            y.extend([miny, maxy])
            hole_set.add(tuple(hole.origin))
        x = np.unique(x)
        y = np.unique(y)

        # LOGGER.info("Expected to be divided {} into {} cuboids.".format(self, (len(x) - 1) * (len(y) - 1)))

        C = list()
        for j in range(len(y) - 1):
            C.append([])
            for i in range(len(x) - 1):
                origin = [x[i], y[j]]
                origin.insert(index, self.origin[index])
                side = [x[i + 1] - x[i], y[j + 1] - y[j]]
                side.insert(index, self.side[index])
                C[j].append(Cuboid(origin=tuple(origin), side=tuple(side), ele_size=self.ele_size).initial())

        # C = [[Cuboid(origin=(x[i], y[j], self.origin[2]),
        #              side=(x[i + 1] - x[i], y[j + 1] - y[j], self.side[2]),
        #              ele_size=self.ele_size).initial() for i in range(len(x) - 1)]
        #      for j in range(len(y) - 1)]

        for j in range(len(y) - 1):
            for i in range(len(x) - 1):
                if tuple(C[j][i].origin) in hole_set:
                    # LOGGER.debug("{} is a hole.".format(C[j][i]))
                    neighbors = [(j + 1, i), (j - 1, i), (j, i - 1), (j, i + 1)]
                    for p, q in neighbors:
                        if -1 < p < len(y) - 1 and -1 < q < len(x) - 1:
                            C[p][q].clarify_bound_with(C[j][i])

        points = [C[j][i].points for j in range(len(y) - 1) for i in range(len(x) - 1)]
        return np.concatenate(points)

    def _bound_points_initial(self):
        if self.holes is not None:
            points = self._divide()
        else:
            points = super(Block, self)._points_initial()
        nx, ny, nz = np.ceil(self.side / self.ele_size)
        dx, dy, dz = [e[0] / e[1] for e in zip(self.side, [nx, ny, nz])]
        temp = (points - self.origin) / np.array([dx, dy, dz])
        left_index = np.argwhere(temp[:, 0] <= 1).flatten()
        right_index = np.argwhere(temp[:, 0] >= nx - 1).flatten()
        front_index = np.argwhere(temp[:, 1] <= 1).flatten()
        back_index = np.argwhere(temp[:, 1] >= ny - 1).flatten()
        down_index = np.argwhere(temp[:, 2] <= 1).flatten()
        up_index = np.argwhere(temp[:, 2] >= nz - 1).flatten()
        index = [left_index, right_index, front_index, back_index, down_index, up_index]
        bound_index = np.unique(np.concatenate(index))
        return points[bound_index]

    def _inner_points_initial(self):
        bound = Cuboid(origin=self.origin + 1.5 * self.ele_size, side=self.side - 3 * self.ele_size)
        holes = None
        if self.holes is not None:
            index = np.argwhere(self.holes[0].side == self.side)[0][0]
            holes = list()
            for hole in self.holes:
                pan = np.zeros(3)
                pan[index] = 1.5 * self.ele_size
                origin = hole.origin + pan
                # origin = hole.origin + np.array([0, 0, 1.5 * self.ele_size])
                side = hole.side - 2 * pan
                # side = hole.side - np.array([0, 0, 3 * self.ele_size])
                holes.append(Cuboid(origin=origin, side=side))
        inner_block = Block(bound=bound, holes=holes, ele_size=2 * self.ele_size)
        if inner_block.holes is not None:
            return inner_block._divide()
        else:
            return super(Block, inner_block)._points_initial()

    def _points_initial(self):
        # if self.holes is None:
        #     return super(Block, self)._points_initial()
        # else:
            return np.concatenate((self._bound_points_initial(), self._inner_points_initial()))

    def pan(self, vector):
        super(Block, self).pan(vector)
        if self.holes is not None:
            for hole in self.holes:
                hole.pan(vector)

    def rotate(self, centroid, direction: str=None):
        # rotate bound and holes
        super(Block, self).rotate(centroid, direction)
        if self.holes is not None:
            for hole in self.holes:
                hole.rotate(centroid, direction)

    def in_hole(self, index):
        def arg(e, arr: list):
            if len(arr) == 0:
                return 0
            elif e < arr[0]:
                return 0
            elif e > arr[-1]:
                return len(arr) - 1
            n = len(arr) // 2
            if arr[n] == e:
                return n
            elif arr[n] < e:
                return n + arg(e, arr[n:])
            else:
                return arg(e, arr[:n])

        p, q = {0, 1, 2} - {index}
        diagonal_point = self.origin + self.side
        x, y = [self.origin[p], diagonal_point[p]], [self.origin[q], diagonal_point[q]]
        hole_set = set()
        if self.holes is not None:
            for hole in self.holes:
                minx, miny = hole.origin[[p, q]]
                maxx, maxy = (hole.origin + hole.side)[[p, q]]
                x.extend([minx, maxx])
                y.extend([miny, maxy])
                hole_set.add((minx, miny))
        x = np.unique(x)
        y = np.unique(y)

        hole_index = set([(i, j) for j in range(len(y) - 1) for i in range(len(x) - 1) if (x[i], y[j]) in hole_set])

        # nz = np.ceil(self.side[2] / self.ele_size)
        # dz = self.side[2] / nz

        def _in_hole(point):
            pointx, pointy = point[[p, q]]
            ix = arg(pointx, x)
            jy = arg(pointy, y)
            return (ix, jy) in hole_index

        return _in_hole

    # @logging("building")
    def build(self):
        LOGGER.info("start building {}...".format(self))
        points = self.points
        n = points.shape[0]
        units = list()

        for plane in self.planes():
            points = np.r_[points, mirror_projection(self.points, plane)]

        # see the {@ref: http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html#scipy.spatial.Voronoi}
        # to learn more about Voronoi and its attributes
        vor = Voronoi(points)

        for i, ith_region in enumerate(vor.point_region[:n]):
            curr_region = vor.regions[ith_region]
            region_map = dict(zip(curr_region, list(range(len(curr_region)))))
            curr_verts = vor.vertices[curr_region]

            ridge_index = np.r_[vor.ridge_points[vor.ridge_points[:, 0] == i],
                                vor.ridge_points[vor.ridge_points[:, 1] == i]]
            ridge = [vor.ridge_dict[tuple(index)]
                     for index in ridge_index if -1 not in vor.ridge_dict[tuple(index)]]

            # mapping ridge
            ridge = [[region_map[e] for e in r] for r in ridge]

            units.append(Unit(verts=curr_verts, ridge=ridge, material=self.material))

        self.units = units
        return self

    def remove_redundancy(self, zlist: set=None, index: int=2):
        pin_down, pin_up = False, False
        if zlist is not None:
            pin_down = self.origin[index] in zlist
            pin_up = self.origin[index] + self.side[index] in zlist
        n = np.ceil(self.side[index] / self.ele_size)
        d = self.side[index] / n
        units = list()
        for i, unit in enumerate(self.units):
            temp = (self.points[i][index] - self.origin[index]) / d
            if not self.in_hole(index)(self.points[i]):
                units.append(unit)
            elif (pin_down and temp <= 1) or (pin_up and temp >= n - 1):
                unit.material = Block.morta
                units.append(unit)
        self.units = units
        return self

    # @logging("plotting")
    def plot(self, figure=None, **kwargs):
        if not hasattr(self, "units"):
            raise TypeError("must build first.")
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        # LOGGER.info("start plotting {}...".format(self))
        figure = Axes3D(plt.figure()) if not figure else figure
        vision_planes = set([str(bound_plane) for bound_plane in self.planes()])
        if self.holes is not None:
            for hole in self.holes:
                vision_planes = vision_planes.union(set([str(hole_plane) for hole_plane in hole.planes()]))

        for unit in self.units:
            for unit_plane in unit.planes():
                if unit.material == Block.brick:
                    if str(unit_plane) in vision_planes:
                        unit_plane.plot(figure, **Block.attribute[unit.material])
                else:
                    unit_plane.plot(figure, **Block.attribute[unit.material])

    def output(self, path):
        pass

    def __repr__(self):
        return "{}:{}->{}".format(self.material, self.origin, self.origin + self.side)


if __name__ == "__main__":
    origin = (0, 0, 0)
    side = (240, 115, 90)
    bound = Cuboid(origin=origin, side=side)
    # model.build()
    # figure = Axes3D(plt.figure())
    # model.plot(figure)
    # minx, maxx = figure.get_xlim()
    # miny, maxy = figure.get_ylim()
    # minz, maxz = figure.get_zlim()
    # max_side = max(side)
    # figure.set_xlim(minx, minx + max_side)
    # figure.set_ylim(miny, miny + max_side)
    # figure.set_zlim(minz, minz + max_side)
    # plt.show()

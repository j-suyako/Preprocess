from common import Unit, Cuboid
from lib.simple3D import *
from util.logging import logging


# TODO: modify the initiation method of points in model
def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.intercept) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


def sandwich(collection, breads: list, axis: int, sort=True):
    """get the cross sections perpendicular to {axis} in {intercepts} in {collection}

    :param collection:
    :param breads:
    :param axis:
    :param sort:
    :return:
    """
    for e in collection:
        if 'range' not in e.__dict__:
            raise TypeError()
    if axis not in [0, 1, 2]:
        raise ValueError()
    if sort:
        collection.sort(key=lambda x: x.range[axis][0])
    breads.sort()
    start_index, end_index = 0, 0
    normal_vector = np.zeros(3)
    normal_vector[axis] = 1
    for bottom, above in breads:
        bottom_plane = InfinitePlane(normal_vector, -bottom)
        above_plane = InfinitePlane(normal_vector, -above)
        res = list()
        if collection[end_index].range[axis][1] < bottom:
            start_index = end_index
        else:
            end_index = start_index
        for unit in collection[start_index:]:
            if unit.range[axis][0] <= bottom < unit.range[axis][1]:
                res.append(unit.cutby(bottom_plane)[1])
            elif unit.range[axis][0] < above <= unit.range[axis][1]:
                res.append(unit.cutby(above_plane)[0])
            elif bottom < unit.range[axis][0] and above > unit.range[axis][1]:
                res.append(unit)
            elif unit.range[axis][0] >= above:
                break
            end_index += 1
        yield res


def split(collection, axis: int, sort=True):
    for e in collection:
        if 'range' not in e.__dict__:
            raise TypeError()
    if axis not in [0, 1, 2]:
        raise ValueError()
    if sort:
        collection.sort()
    res = list()
    res.append(0)
    n = len(collection)
    for i in range(n - 1):
        if collection[i].range[axis][0] < collection[i + 1].range[axis][0]:
            res.append(i + 1)
    res.append(n)
    return res


def get_range(collection, slices, axis: int):
    return [(collection[s].range[axis][0], collection[s].range[axis][1]) for s in slices[:-1]]


class Model(object):
    def __init__(self, bound: Cuboid, ele_size=5.0, holes=None):
        """

        :param bound:
        :param ele_size:
        :param holes:
        """
        self.bound = bound
        self.ele_size = ele_size
        self.holes = holes
        self.points = self._points_initial()
        self.units = list()

    @logging("points initialization")
    def _points_initial(self) -> np.ndarray:
        """points initialization

        :return:
        """
        # TODO: judge if points in bound after initialization(irregular bound)
        n_points = int(0.68 * self.bound.volume() / self.ele_size ** 3)

        zoom_matrix = np.zeros((3, 3), dtype="float64")
        for i in range(3):
            zoom_matrix[i][i] = self.bound.range[i][1] - self.bound.range[i][0]

        pan_matrix = np.array([self.bound.range[0][0], self.bound.range[1][0], self.bound.range[2][0]], dtype="float64")
        return np.dot(np.random.rand(n_points, 3), zoom_matrix) + pan_matrix

    @logging("quasihomogenization")
    def _quasi_homogenization(self, iter_num=6):
        """first homogenization the distribution of points

        :param iter_num:
        :return:
        """
        n = self.points.shape[0]

        LOGGER.info("start homogenization...")
        for i in range(iter_num):
            LOGGER.debug("{}th homogenization".format(i + 1))
            points = self.points
            for plane in self.bound.planes():
                points = np.r_[points, mirror_projection(self.points, plane)]

            # see the {@ref: http://scipy.github.io/devdocs/generated/scipy.spatial.Voronoi.html#scipy.spatial.Voronoi}
            # to learn more about Voronoi and its attributes
            ith_vor = Voronoi(points)
            points = np.array([])

            for j, jth_region in enumerate(ith_vor.point_region[:n]):
                curr_region = ith_vor.regions[jth_region]
                region_map = dict(zip(curr_region, list(range(len(curr_region)))))
                curr_points = ith_vor.vertices[curr_region]

                ridge_index = np.r_[ith_vor.ridge_points[ith_vor.ridge_points[:, 0] == j],
                                    ith_vor.ridge_points[ith_vor.ridge_points[:, 1] == j]]
                ridge = [ith_vor.ridge_dict[tuple(index)] for index in ridge_index if
                         -1 not in ith_vor.ridge_dict[tuple(index)]]

                # mapping ridge
                ridge = [[region_map[e] for e in r] for r in ridge]

                jth_unit = Unit(index=j, points=curr_points, ridge=ridge)
                points = np.r_[points, jth_unit.centroid()]
                if i == iter_num - 1:
                    self.units.append(jth_unit)

            self.points = points.reshape((-1, 3))

    @logging("build")
    def build(self, quasi_homo_iter=6):
        self._quasi_homogenization(quasi_homo_iter)
        if self.holes is not None:
            LOGGER.info("opening hole...")
            self.assign(self.holes, Unit.hole)
        self.units = [e for e in self.units if e.material != Unit.hole]

    @logging("assigning")
    def assign(self, regions: [list, str], material: str):
        if isinstance(regions, list):
            LOGGER.info("assigning {} by {}".format(regions, material))
            in_region_units = list()
            bread_slice = split(collection=regions, axis=2)
            bread_range = get_range(collection=regions, slices=bread_slice, axis=2)
            for i, meat in enumerate(sandwich(self.units, bread_range, 2)):  # meat layer
                LOGGER.debug("get meat layer between {} with {} in z axis".format(bread_range[i][0], bread_range[i][1]))
                meat_region = regions[bread_slice[i]:bread_slice[i + 1]]
                meat_slice = split(collection=meat_region, axis=0, sort=False)
                meat_range = get_range(collection=meat_region, slices=meat_slice, axis=0)
                for j, butter in enumerate(sandwich(meat, meat_range, 0)):  # butter layer
                    LOGGER.debug("get butter layer between {} with {} in x axis".format(meat_range[j][0], meat_range[j][1]))
                    butter_region = meat_region[meat_slice[j]:meat_slice[j + 1]]
                    butter_slice = split(collection=butter_region, axis=1, sort=False)
                    butter_range = get_range(collection=butter_region, slices=butter_slice, axis=1)
                    for k, lettuce in enumerate(sandwich(butter, butter_range, 1)):
                        LOGGER.debug("get lettuce between {} with {} in y axis".format(butter_range[k][0], butter_range[k][1]))
                        in_region_units.append([e.index for e in lettuce])
            self.units.sort(key=lambda x: x.index)
            for i, region in enumerate(regions):
                if not isinstance(region, Cuboid):
                    raise TypeError()
                LOGGER.debug("{}th region cutting...".format(i))
                for j in in_region_units[i]:
                    unit = self.units[j]
                    try:
                        in_poly, out_poly = region.cut(unit)
                        if in_poly is None:
                            continue
                        elif out_poly is None:
                            unit.material = material
                        else:
                            self.units[j] = Unit(index=unit.index, poly=out_poly, material=unit.material)
                            self.units.append(Unit(index=len(self.units), poly=in_poly, material=material))
                    except Exception:
                        raise ValueError("\nthe region is:\n{}\nthe point of unit is:\n{}".format(region, unit.points))
        elif isinstance(regions, str):
            if regions == "remain":
                for unit in self.units:
                    if unit.material is None:
                        unit.material = material

    @logging("plotting")
    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        LOGGER.info("start plotting...")
        figure = Axes3D(plt.figure()) if not figure else figure
        vision_planes = set([str(bound_plane) for bound_plane in self.bound.planes()])
        if self.holes is not None:
            for hole in self.holes:
                vision_planes.union(set([str(hole_plane) for hole_plane in hole.planes()]))
        for unit in self.units:
            if unit.material != Unit.hole:
                for unit_plane in unit.planes():
                    if str(unit_plane) in vision_planes:
                        if unit.material == Unit.brick:
                            attribute = {"alpha": 0.5, "color": (0.42, 0.1, 0.05)}
                        elif unit.material == Unit.mortar:
                            attribute = {"alpha": 0.5, "color": (0.75, 0.75, 0.75)}
                        else:
                            attribute = {"alpha": 0.5, "color": (0.42, 0.1, 0.05)}
                        unit_plane.plot(figure, **attribute)

    def output(self, path: str):
        pass


if __name__ == "__main__":
    origin = (0, 0, 0)
    side = (240, 115, 90)
    bound = Cuboid(origin=origin, side=side)
    holes = list()
    hole1 = Cuboid(origin=(15, 15, 0), side=(35, 35, 90))
    # hole2 = Cuboid(origin=(15, 65, 0), side=(35, 35, 90))
    # hole4 = Cuboid(origin=(65, 65, 0), side=(47.5, 35, 90))
    # hole3 = Cuboid(origin=(65, 15, 0), side=(47.5, 35, 90))
    # hole5 = Cuboid(origin=(127.5, 15, 0), side=(47.5, 35, 90))
    # hole6 = Cuboid(origin=(127.5, 65, 0), side=(47.5, 35, 90))
    # hole7 = Cuboid(origin=(190, 15, 0), side=(35, 35, 90))
    # hole8 = Cuboid(origin=(190, 65, 0), side=(35, 35, 90))
    # hole2 = Cuboid(origin=(15, 140, 0), side=(35, 35, 90))
    # hole3 = Cuboid(origin=(265, 15, 0), side=(35, 35, 90))
    holes.append(hole1)
    # holes.append(hole2)
    # holes.append(hole3)
    # holes.append(hole2)
    # holes.append(hole3)
    # holes.append(hole4)
    # holes.append(hole5)
    # holes.append(hole6)
    # holes.append(hole7)
    # holes.append(hole8)
    model = Model(bound=bound, ele_size=20, holes=holes)
    model.build()
    # brick_bound1 = Cuboid(origin=(0, 0, 0), side=(240, 115, 90))
    # brick_bound2 = Cuboid(origin=(0, 125, 0), side=(240, 115, 90))
    # brick_bound3 = Cuboid(origin=(250, 0, 0), side=(115, 240, 90))
    # model.assign([brick_bound1, brick_bound2, brick_bound3], Unit.brick)
    # model.assign("remain", Unit.mortar)
    figure = Axes3D(plt.figure())
    model.plot(figure)
    minx, maxx = figure.get_xlim()
    miny, maxy = figure.get_ylim()
    minz, maxz = figure.get_zlim()
    max_side = max(side)
    figure.set_xlim(minx, minx + max_side)
    figure.set_ylim(miny, miny + max_side)
    figure.set_zlim(minz, minz + max_side)
    plt.show()

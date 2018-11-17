from lib.simple3D import *
from common import Unit, Cuboid


# TODO: modify the initiation method of points in model
def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.intercept) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


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

    def _points_initial(self) -> np.ndarray:
        """
        points initialization
        :return:
        """
        # TODO: judge if points in bound after initialization(irregular bound)
        n_points = int(0.68 * self.bound.volume() / self.ele_size ** 3)

        zoom_matrix = np.array([[0] * 3 for _ in range(3)], dtype="float64")
        for i in range(3):
            zoom_matrix[i][i] = self.bound.range[i][1] - self.bound.range[i][0]

        pan_matrix = np.array([self.bound.range[0][0], self.bound.range[1][0], self.bound.range[2][0]], dtype="float64")
        return np.dot(np.random.rand(n_points, 3), zoom_matrix) + pan_matrix

    def _quasi_homogenization(self, iter_num=6):
        """
        first homogenization the distribution of points
        :param iter_num:
        :return:
        """
        n = self.points.shape[0]

        for i in range(iter_num):
            points = self.points
            for r in self.bound.ridge:
                plane = Plane(self.bound.points[r])
                points = np.r_[points, mirror_projection(points, plane)]

            ith_vor = Voronoi(points)
            points = np.array([])

            for j, jth_region in enumerate(ith_vor.point_region[:n]):
                curr_region = ith_vor.regions[jth_region]
                region_map = dict(zip(curr_region, list(range(len(curr_region)))))
                curr_points = ith_vor.vertices[curr_region]

                # if np.min(curr_points) > -1e-3 and np.max(curr_points) < 1 + 1e-3:
                ridge_index = np.r_[ith_vor.ridge_points[ith_vor.ridge_points[:, 0] == j],
                                    ith_vor.ridge_points[ith_vor.ridge_points[:, 1] == j]]
                ridge = [ith_vor.ridge_dict[tuple(index)] for index in ridge_index if
                         -1 not in ith_vor.ridge_dict[tuple(index)]]

                # mapping ridge
                for _, r in enumerate(ridge):
                    ridge[_] = [region_map[e] for e in r]

                jth_unit = Unit(points=curr_points, ridge=ridge)
                points = np.r_[points, jth_unit.centroid()]
                if i == iter_num - 1:
                    self.units.append(jth_unit)

            self.points = points.reshape((-1, 3))

    def build(self, quasi_homo_iter=6):
        self._quasi_homogenization(quasi_homo_iter)
        if self.holes is not None:
            self.assign(self.holes, Unit.hole)

    def assign(self, regions: [list, str], material: str):
        if isinstance(regions, list):
            for region in regions:
                units = list()
                for unit in self.units:
                    in_region, out_region = region.cut(unit)
                    if in_region is not None and material != Unit.hole:
                        units.append(Unit(in_region, material))
                    if out_region is not None:
                        units.append(Unit(out_region))
                self.units = units
        elif isinstance(regions, str):
            if regions == "remain":
                for unit in self.units:
                    if unit.material is None:
                        unit.material = material

    def plot(self, figure=None):
        if figure and not isinstance(figure, Axes3D):
            raise TypeError("figure argument must be a 3d layer")
        figure = Axes3D(plt.figure()) if not figure else figure
        vision_planes = set([str(bound_plane) for bound_plane in self.bound.planes()])
        if self.holes is not None:
            for hole in self.holes:
                vision_planes.union(set([str(hole_plane) for hole_plane in hole.planes()]))
        for unit in self.units:
            for unit_plane in unit.planes():
                # if str(unit_plane) in vision_planes:
                    unit_plane.plot(figure)

    def output(self, path: str):
        pass


if __name__ == "__main__":
    bound = Cuboid(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]))
    holes = list()
    hole = Cuboid(np.array([[0.2, 0.2, 0], [0.8, 0.2, 0], [0.8, 0.8, 0], [0.2, 0.8, 0],
                            [0.8, 0.8, 1], [0.8, 0.2, 1], [0.8, 0.8, 1], [0.2, 0.8, 1]]))
    holes.append(hole)
    model = Model(bound=bound, ele_size=0.2, holes=holes)
    model.build()
    model.plot()
    plt.show()

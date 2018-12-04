from lib.simple3D import *
from util.math import continuous_block


class Cuboid(Polyhedron):
    def __init__(self, origin, side):
        x, y, z = origin
        a, b, c = side
        points = np.array([[px, py, pz] for px in [x, x + a] for py in [y, y + b] for pz in [z, z + c]])
        ridge = [[0, 1, 3, 2], [0, 1, 5, 4], [0, 2, 6, 4], [1, 3, 7, 5], [2, 3, 7, 6], [4, 5, 7, 6]]
        super(Cuboid, self).__init__(points, ridge)
        self.origin = origin
        self.end = (x + a, y + b, z + c)

    def volume(self):
        return np.prod([self.range[i][1] - self.range[i][0] for i in range(3)])

    def centroid(self):
        return (self.range[:, 0] + self.range[:, 1]) / 2

    def contains(self, points):
        if len(points.shape) == 1:
            around_point = np.around(points, decimals=5)
            return np.all([self.range[i][0] <= around_point[i] <= self.range[i][1] for i in range(3)])
        else:
            return [self.contains(point) for point in points]

    def cut(self, poly: Polyhedron):
        is_contain = np.array(self.contains(poly.points))

        if not np.all(is_contain) and not np.all(~is_contain):
            intersect_planes = list()
            for plane in self.planes():
                dis_arr = plane.dis_to(poly.points)
                if not np.all(dis_arr >= 0) and not np.all(dis_arr <= 0):
                    points = poly.intersect(InfinitePlane(normal_vector=plane.normal_vector, intercept=plane.intercept))
                    if not np.all(~np.array(plane.contains(points))):
                        intersect_planes.append(plane)

            try:
                if len(intersect_planes) == 1:
                    return self._cut_1_plane(poly, intersect_planes[0])
                elif len(intersect_planes) == 2:
                    return self._cut_2_plane(poly, intersect_planes, is_contain)
                elif len(intersect_planes) == 3:
                    return self._cut_3_plane(poly, intersect_planes, is_contain)
                elif len(intersect_planes) > 3:
                    raise ValueError()
            except Exception:
                raise ValueError("in points are:\n{}\n\nout points are:\n{}".format(poly.points[is_contain],
                                                                                    poly.points[~is_contain]))

        elif np.all(is_contain):
            return poly, None

        else:  # len(out_points) > 0
            return None, poly

    def _cut_1_plane(self, poly: Polyhedron, plane: Plane):
        poly1, poly2 = poly.cutby(plane)
        if self.contains(poly1.points[0]):
            return poly1, poly2
        else:
            return poly2, poly1

    def _cut_2_plane(self, poly: Polyhedron, planes, is_contain):
        plane0, plane1 = planes[0], planes[1]
        line = set([tuple(e) for e in plane0.points]) & set([tuple(e) for e in plane1.points])
        if len(line) != 2:
            raise ValueError("ele size seems too large")

        intersect_points, keys = list(), list()
        for plane in planes:
            points, key = poly.intersect(plane, return_key=True)
            intersect_points.append(points)
            keys.extend(key)
        for line in [line]:
            points, key = poly.intersect(Line(line.pop(), line.pop()), return_key=True)
            if len(key) == 0:
                return None, None
            intersect_points.append(points)
            keys.extend(key)

        n = poly.points.shape[0]
        mapping = dict()
        for i, key in enumerate(keys):
            mapping.setdefault(key, list())
            mapping[key].append(i + n)
        # mapping = dict(zip(keys, list(range(n, n + len(keys)))))

        ridge = [[], []]
        for r in poly.ridge:
            # in_flag = list()
            new_r = list()
            for i in range(len(r)):
                start, end = i, i + 1 if i < len(r) - 1 else 0
                minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                new_r.append(r[start])
                if (minr, maxr) in mapping:
                    # in_flag.append(start)
                    new_r.extend(mapping[(minr, maxr)])
            in_flag = list(np.where([1 if e >= n or is_contain[e] else 0 for e in new_r])[0])
            if np.size(in_flag) == len(new_r):
                ridge[0].append(new_r)
            elif np.size(in_flag) > 2:
                ridge[0].append(list(np.array(new_r)[in_flag]))
                block = continuous_block(list(in_flag))
                for i in range(len(block) - 1):
                    ridge[1].append(new_r[block[i][-1]:block[i + 1][0] + 1])
                if len(block) > 0 and (block[0][0] != 0 or block[-1][-1] != len(new_r) - 1):
                    ridge[1].append(new_r[block[-1][-1]:] + new_r[:block[0][0] + 1])
            else:
                ridge[1].append(new_r)
            # if len(in_flag) < 2:
            #     ridge[0 if self.contains(poly.points[r[0]]) else 1].append(r)
            # else:
            #     r = tuple(r)
            #     if r in mapping:
            #         r0, r1 = new_r[:in_flag[0] + 2] + [mapping[r]] + new_r[in_flag[1] + 2:], new_r[in_flag[0] + 1:in_flag[1] + 3] + [mapping[r]]
            #     else:
            #         r0, r1 = new_r[:in_flag[0] + 2] + new_r[in_flag[1] + 2:], new_r[in_flag[0] + 1:in_flag[1] + 3]
            #     if not self.contains(poly.points[r[0]]):
            #         r0, r1 = r1, r0
            #     ridge[0].append(r0)
            #     ridge[1].append(r1)

        r0 = FinitePlane(np.r_[intersect_points[0], intersect_points[2]]).ridge
        for i, e in enumerate(r0):
            remap = n
            if e >= intersect_points[0].shape[0]:
                remap += intersect_points[1].shape[0]
            r0[i] = e + remap

        r1 = FinitePlane(np.r_[intersect_points[1], intersect_points[2]]).ridge
        for i, e in enumerate(r1):
            remap = n + intersect_points[0].shape[0]
            r1[i] = e + remap

        for r in [r0, r1]:
            ridge[0].append(list(r))
            ridge[1].append(list(r))

        intersect_points = np.concatenate([e for e in intersect_points])
        in_points = np.r_[poly.points[is_contain], intersect_points]
        out_points = np.r_[poly.points[~is_contain], intersect_points]

        in_map = list(np.arange(n)[is_contain])
        out_map = list(np.arange(n)[~is_contain])
        in_map.extend(list(range(n, n + len(keys))))
        out_map.extend(list(range(n, n + len(keys))))

        return Polyhedron(in_points, ridge=ridge_renum(in_map, ridge[0])), \
               Polyhedron(out_points, ridge=ridge_renum(out_map, ridge[1]))

    def _cut_3_plane(self, poly: Polyhedron, planes, is_contain):
        plane0, plane1, plane2 = planes[0], planes[1], planes[2]
        line0 = set([tuple(e) for e in plane0.points]) & set([tuple(e) for e in plane1.points])
        line1 = set([tuple(e) for e in plane0.points]) & set([tuple(e) for e in plane2.points])
        line2 = set([tuple(e) for e in plane1.points]) & set([tuple(e) for e in plane2.points])
        corner_point = line0 & line1 & line2
        if len(corner_point) != 1:
            raise ValueError()

        intersect_points, keys = list(), list()
        for plane in planes:
            points, key = poly.intersect(plane, return_key=True)
            intersect_points.append(points)
            keys.extend(key)
        for line in [line0, line1, line2]:
            points, key = poly.intersect(Line(line.pop(), line.pop()), return_key=True)
            intersect_points.append(points)
            keys.extend(key)
        intersect_points.append(np.array(list(corner_point)))

        n = poly.points.shape[0]
        mapping = dict(zip(keys, list(range(n, n + len(keys)))))

        ridge = [[], []]
        for r in poly.ridge:
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
                ridge[0 if self.contains(poly.points[r[0]]) else 1].append(r)
            else:
                r = tuple(r)
                if r in mapping:
                    r0, r1 = new_r[:flag[0] + 2] + [mapping[r]] + new_r[flag[1] + 2:], new_r[flag[0] + 1:flag[1] + 3] + [mapping[r]]
                else:
                    r0, r1 = new_r[:flag[0] + 2] + new_r[flag[1] + 2:], new_r[flag[0] + 1:flag[1] + 3]
                if not self.contains(poly.points[r[0]]):
                    r0, r1 = r1, r0
                ridge[0].append(r0)
                ridge[1].append(r1)

        r0 = FinitePlane(np.r_[intersect_points[0], intersect_points[3], intersect_points[4], intersect_points[-1]]).ridge
        for i, e in enumerate(r0):
            remap = n
            if e >= intersect_points[0].shape[0]:
                remap += intersect_points[1].shape[0] + intersect_points[2].shape[0]
            if e >= intersect_points[0].shape[0] + intersect_points[3].shape[0] + intersect_points[4].shape[0]:
                remap += intersect_points[5].shape[0]
            r0[i] = e + remap

        r1 = FinitePlane(np.r_[intersect_points[1], intersect_points[3], intersect_points[5], intersect_points[-1]]).ridge
        for i, e in enumerate(r1):
            remap = n + intersect_points[0].shape[0]
            if e >= intersect_points[1].shape[0]:
                remap += intersect_points[2].shape[0]
            if e >= intersect_points[1].shape[0] + intersect_points[3].shape[0]:
                remap += intersect_points[4].shape[0]
            r1[i] = e + remap

        r2 = FinitePlane(np.r_[intersect_points[2], intersect_points[4], intersect_points[5], intersect_points[-1]]).ridge
        for i, e in enumerate(r2):
            remap = n + intersect_points[0].shape[0] + intersect_points[1].shape[0]
            if e >= intersect_points[2].shape[0]:
                remap += intersect_points[3].shape[0]
            r2[i] = e + remap

        for r in [r0, r1, r2]:
            ridge[0].append(list(r))
            ridge[1].append(list(r))

        intersect_points = np.concatenate([e for e in intersect_points])
        in_points = np.r_[poly.points[is_contain], intersect_points]
        out_points = np.r_[poly.points[~is_contain], intersect_points]

        in_map = list(np.arange(n)[is_contain])
        out_map = list(np.arange(n)[~is_contain])
        in_map.extend(list(range(n, n + len(keys) + 1)))
        out_map.extend(list(range(n, n + len(keys) + 1)))

        return Polyhedron(in_points, ridge=ridge_renum(in_map, ridge[0])), \
               Polyhedron(out_points, ridge=ridge_renum(out_map, ridge[1]))

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
        return "{} -> {}".format(self.origin, self.end)


if __name__ == "__main__":
    hole = Cuboid(origin=(680, 0, 400), side=(640, 240, 500))
    # holes = [hole]
    # print(holes)
    points = np.array([[816.18987, 92.98193, 422.85506],
                       [817.66813, 0.0, 427.84348],
                       [689.77081, 0.0, 460.98953],
                       [833.21866, 65.40623, 329.62482],
                       [687.38827, 72.51114, 457.41806],
                       [651.68991, 92.29704, 413.50438],
                       [818.86619, 104.04854, 405.9786],
                       [824.14898, 97.92449, 376.70507],
                       [702.57177, 133.83459, 398.80186],
                       [711.03524, 116.85628, 307.43827],
                       [783.56547, 117.29774, 429.90533],
                       [752.52741, 126.76453, 437.4023],
                       [790.46602, 93.79913, 295.50935],
                       [810.71212, 91.05519, 311.90677],
                       [624.31278, 58.89258, 380.26408],
                       [632.07646, 47.32298, 311.36241],
                       [623.72706, -0.0, 380.07645],
                       [632.05283, 0.0, 307.30435],
                       [784.14343, 0.0, 288.252],
                       [834.89333, 0.0, 329.50986]])
    ridge = [[0, 6, 10], [1, 19, 3, 7, 6, 0], [2, 16, 14, 5, 4], [16, 17, 15, 14], [17, 18, 12, 9, 15], [18, 12, 13, 3, 19], [3, 7, 13], [4, 11, 8, 5], [5, 14, 15, 9, 8], [6, 7, 13, 12, 9, 8, 11, 10], [0, 1, 2, 4, 11, 10], [2, 1, 19, 18, 17, 16]]
    poly = Polyhedron(points=points, ridge=ridge)
    # poly2 = Polyhedron(points=points, ridge=ridge)
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    poly.plot(ax1, alpha=0.2)
    hole.plot(ax1, alpha=0.2)
    in_hole, out_hole = hole.cut(poly)
    # in_hole.plot(ax2, alpha=0.2)
    out_hole.plot(ax2, alpha=0.2)
    plt.show()

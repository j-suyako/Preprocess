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

        # if not np.all(is_contain) and not np.all(~is_contain):
        intersect_planes = list()
        for plane in self.planes():
            dis_arr = plane.dis_to(poly.points)
            if not np.all(dis_arr >= 0) and not np.all(dis_arr <= 0):
                points = poly.intersect(InfinitePlane(normal_vector=plane.normal_vector, intercept=plane.intercept))
                if not np.all(~np.array(plane.contains(points))):
                    intersect_planes.append(plane)

        try:
            if len(intersect_planes) == 0:
                if np.all(is_contain):
                    return poly, None
                else:
                    return None, poly
            elif len(intersect_planes) == 1:
                return self._cut_1_plane(poly, intersect_planes[0])
            elif len(intersect_planes) == 2:
                return self._cut_2_plane(poly, intersect_planes, is_contain)
            elif len(intersect_planes) == 3:
                # return None, None
                return self._cut_3_plane(poly, intersect_planes, is_contain)
            elif len(intersect_planes) > 3:
                raise ValueError()
        except Exception:
            return None, None
            # raise ValueError("in points are:\n{}\n\nout points are:\n{}".format(poly.points[is_contain],
            #                                                                     poly.points[~is_contain]))

        # elif np.all(is_contain):
        #     return poly, None
        #
        # else:  # len(out_points) > 0
        #     return None, poly

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
            points, key = poly.intersect(Segment(line.pop(), line.pop()), return_key=True)
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

        all_points = np.concatenate((poly.points, np.concatenate([e for e in intersect_points])))
        ridge = [[], []]
        for r in poly.ridge:
            if tuple(r) in mapping:
                in_flag = list()
                new_r = list()
                for i in range(len(r)):
                    start, end = i, i + 1 if i < len(r) - 1 else 0
                    minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                    new_r.append(r[start])
                    if (minr, maxr) in mapping:
                        in_flag.append(start)
                        t = mapping[(minr, maxr)]
                        if len(t) == 2:
                            vec0 = all_points[t[0]] - all_points[r[start]]
                            vec1 = all_points[t[1]] - all_points[r[start]]
                            length0 = np.dot(vec0, vec0)
                            length1 = np.dot(vec1, vec1)
                            if length0 > length1:
                                new_r.extend([t[1], t[0]])
                            else:
                                new_r.extend(t)
                        else:
                            new_r.extend(t)
                r = tuple(r)
                if len(in_flag) < 2:
                    r0, r1 = new_r[:in_flag[0] + 2] + mapping[r] + new_r[in_flag[0] + 2:], new_r[in_flag[0] + 1:in_flag[0] + 3] + mapping[r]
                else:
                    r0, r1 = new_r[:in_flag[0] + 2] + mapping[r] + new_r[in_flag[1] + 2:], new_r[in_flag[0] + 1:in_flag[1] + 3] + mapping[r]
                if not self.contains(poly.points[r[0]]):
                    r0, r1 = r1, r0
                ridge[0].append(r0)
                ridge[1].append(r1)
            else:
                new_r = list()
                for i in range(len(r)):
                    start, end = i, i + 1 if i < len(r) - 1 else 0
                    minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                    new_r.append(r[start])
                    if (minr, maxr) in mapping:
                        # in_flag.append(start)
                        t = mapping[(minr, maxr)]
                        if len(t) == 2:
                            vec0 = all_points[t[0]] - all_points[r[start]]
                            vec1 = all_points[t[1]] - all_points[r[start]]
                            length0 = np.dot(vec0, vec0)
                            length1 = np.dot(vec1, vec1)
                            if length0 > length1:
                                new_r.extend([t[1], t[0]])
                            else:
                                new_r.extend(t)
                        else:
                            new_r.extend(t)
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
        if not np.any(self.contains(poly.points)):
            return None, None
        # if not poly.contains(corner_point):
        #     return None, None
        intersect_points, keys = list(), list()
        for plane in planes:
            points, key = poly.intersect(plane, return_key=True)
            intersect_points.append(points)
            keys.extend(key)
        for line in [line0, line1, line2]:
            points, key = poly.intersect(Segment(line.pop(), line.pop()), return_key=True)
            intersect_points.append(points)
            keys.extend(key)
        intersect_points.append(np.array(list(corner_point)))

        n = poly.points.shape[0]
        mapping = dict()
        for i, key in enumerate(keys):
            mapping.setdefault(key, list())
            mapping[key].append(i + n)

        all_points = np.concatenate((poly.points, np.concatenate([e for e in intersect_points])))
        ridge = [[], []]
        for r in poly.ridge:
            in_flag = list()
            new_r = list()
            for i in range(len(r)):
                start, end = i, i + 1 if i < len(r) - 1 else 0
                minr, maxr = (r[start], r[end]) if r[start] < r[end] else (r[end], r[start])
                new_r.append(r[start])
                if (minr, maxr) in mapping:
                    in_flag.append(start)
                    t = mapping[(minr, maxr)]
                    if len(t) == 2:
                        vec0 = all_points[t[0]] - all_points[r[start]]
                        vec1 = all_points[t[1]] - all_points[r[start]]
                        length0 = np.dot(vec0, vec0)
                        length1 = np.dot(vec1, vec1)
                        if length0 > length1:
                            new_r.extend([t[1], t[0]])
                        else:
                            new_r.extend(t)
                    else:
                        new_r.extend(t)
            if len(in_flag) < 2:
                ridge[0 if self.contains(poly.points[r[0]]) else 1].append(r)
            else:
                r = tuple(r)
                if r in mapping:
                    r0, r1 = new_r[:in_flag[0] + 2] + mapping[r] + new_r[in_flag[1] + 2:], new_r[in_flag[0] + 1:in_flag[1] + 3] + mapping[r]
                else:
                    r0, r1 = new_r[:in_flag[0] + 2] + new_r[in_flag[1] + 2:], new_r[in_flag[0] + 1:in_flag[1] + 3]
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
    hole = Cuboid(origin=(65, 65, 0), side=(47.5, 35, 90))
    # holes = [hole]
    # print(holes)
    # points = np.array([[39.28709, 45.45487, 78.01736],
    #                    [37.44342, 43.29013, 79.40019],
    #                    [35.50332, 50.5173, 100.0],
    #                    [64.15525, 28.6708, 88.87884],
    #                    [64.22604, 28.35987, 100.0],
    #                    [67.06721, 31.93641, 90.51979],
    #                    [66.96031, 31.48049, 100.0],
    #                    [34.24821, 48.57623, 100.0],
    #                    [47.84081, 26.01593, 100.0],
    #                    [43.10646, 55.35273, 86.05597],
    #                    [43.31296, 53.04817, 80.77103],
    #                    [52.01929, 47.86635, 75.40178],
    #                    [61.38306, 47.86672, 82.80327],
    #                    [51.38823, 58.48231, 100.0],
    #                    [57.7315, 56.36167, 100.0],
    #                    [47.65265, 26.3383, 87.9009],
    #                    [58.53266, 27.99081, 84.57208],
    #                    [40.83577, 37.66126, 77.51998],
    #                    [46.71891, 40.89292, 73.36141],
    #                    [47.44043, 40.96445, 73.17905]])
    points = np.array([[72.39822, 96.51942, 85.35513],
                       [68.00106, 88.59854, 85.38777],
                       [52.86873, 115.0, 82.04891],
                       [46.02393, 115.0, 100.0],
                       [67.02831, 115.0, 83.58426],
                       [72.26459, 115.0, 100.0],
                       [75.51453, 102.50047, 100.0],
                       [66.60147, 86.44316, 100.0],
                       [48.88944, 90.56151, 83.18921],
                       [52.21523, 87.38798, 83.75394],
                       [53.37971, 85.23931, 100.0],
                       [43.11217, 94.9241, 100.0]])
    # ridge = [[4, 8, 15, 16, 3], [4, 6, 5, 3], [0, 18, 19, 11, 10], [3, 5, 12, 11, 19, 16], [15, 17, 18, 19, 16], [13, 14, 12, 11, 10, 9], [2, 13, 9], [6, 14, 12, 5], [2, 7, 1, 0, 10, 9], [4, 8, 7, 2, 13, 14, 6], [7, 1, 17, 15, 8], [0, 18, 17, 1]]
    ridge = [[4, 2, 8, 9, 1, 0], [5, 6, 7, 10, 11, 3], [4, 5, 6, 0], [6, 7, 1, 0], [2, 4, 5, 3], [11, 10, 9, 8], [7, 10, 9, 1], [2, 3, 11, 8]]
    poly = Polyhedron(points=points, ridge=ridge)
    # poly2 = Polyhedron(points=points, ridge=ridge)
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    poly.plot(ax1, alpha=0.2)
    hole.plot(ax1, alpha=0.2)
    in_hole, out_hole = hole.cut(poly)
    # if in_hole is not None:
    #     in_hole.plot(ax2, alpha=0.2)
    # if out_hole is not None:
    #     out_hole.plot(ax2, alpha=0.2)
    plt.show()

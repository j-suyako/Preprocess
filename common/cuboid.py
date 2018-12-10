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
                return self._cut_3_plane(poly, intersect_planes, is_contain)
            elif len(intersect_planes) > 3:
                raise ValueError()
        except Exception:
            raise ValueError("in points are:\n{}\n\nout points are:\n{}".format(poly.points[is_contain],
                                                                                poly.points[~is_contain]))

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
    hole = Cuboid(origin=(65, 15, 0), side=(47.5, 35, 90))
    # holes = [hole]
    # print(holes)
    points = np.array([[113.69486, 42.95346, 41.85709],
                       [102.71727, 51.11754, 33.26388],
                       [95.83617, 46.42959, 48.28318],
                       [101.92943, 49.62136, 43.9067],
                       [112.07645, 28.00554, 47.34064],
                       [107.08356, 30.37575, 52.67188],
                       [105.62027, 41.31443, 51.78646],
                       [108.53197, 45.00889, 47.46435],
                       [88.54965, 45.60749, 25.66459],
                       [83.29382, 42.62559, 31.56749],
                       [88.90811, 43.57757, 46.03256],
                       [93.25398, 47.67067, 26.01554],
                       [101.11412, 24.30969, 50.86089],
                       [84.32185, 25.31197, 45.1246],
                       [97.60511, 27.63983, 21.19723],
                       [106.23854, 27.32506, 20.87991],
                       [104.79604, 19.7607, 44.35336],
                       [110.85408, 24.40654, 44.80228],
                       [100.35262, 14.8276, 33.51131],
                       [106.43435, 18.79946, 29.19216],
                       [90.73672, 16.6361, 32.4409],
                       [89.74139, 17.8269, 31.34994],
                       [83.40851, 24.4681, 43.20119],
                       [81.66279, 26.76585, 37.08548],
                       [112.3364, 45.30957, 33.51585],
                       [108.02586, 48.43891, 30.55484],
                       [107.47592, 32.55523, 22.07105],
                       [104.34716, 45.34242, 25.15684]])
    ridge = [[0, 24, 25, 1, 3, 7], [0, 4, 17, 19, 15, 26, 24], [1, 25, 27, 11], [4, 17, 16, 12, 5], [8, 14, 15, 26, 27, 11], [12, 13, 22, 20, 18, 16], [8, 14, 21, 23, 9], [16, 18, 19, 17], [20, 22, 23, 21], [14, 15, 19, 18, 20, 21], [0, 4, 5, 6, 7], [2, 6, 5, 12, 13, 10], [1, 3, 2, 10, 9, 8, 11], [9, 23, 22, 13, 10], [2, 6, 7, 3], [24, 26, 27, 25]]
    poly = Polyhedron(points=points, ridge=ridge)
    # poly2 = Polyhedron(points=points, ridge=ridge)
    ax1 = plt.subplot2grid((1, 2), (0, 0), projection="3d")
    ax2 = plt.subplot2grid((1, 2), (0, 1), projection="3d")
    poly.plot(ax1, alpha=0.2)
    hole.plot(ax1, alpha=0.2)
    in_hole, out_hole = hole.cut(poly)
    if in_hole is not None:
        in_hole.plot(ax2, alpha=0.2)
    if out_hole is not None:
        out_hole.plot(ax2, alpha=0.2)
    plt.show()

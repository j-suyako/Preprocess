from Lib.GeometryLib import *


class EhVornoi(Voronoi):
    pass


def mirror_projection(points, face: Plane):
    if len(points.shape) == 1:
        t = -(np.dot(face.normal_vector, points) + face.z_) / (np.dot(face.normal_vector, face.normal_vector))
        return list(points + 2 * t * face.normal_vector)
    else:
        return np.array([mirror_projection(point, face) for point in points])


if __name__ == '__main__':
    # np.random.seed(5)
    # points = np.random.rand(10, 3)
    # polys = [0] * points.shape[0]
    # vors = Voronoi(points)
    # # for pair in vors.ridge_dict:
    # #     points_in_plane = vors.ridge_dict[pair]
    # #     target_point_index = points_in_plane[0] if points_in_plane[0] != -1 else points_in_plane[1]
    # #     this, that = pair
    # #     normal_vector = norm(points[this, :] - points[that, :])
    # #     z_ = -np.dot(normal_vector, points[target_point_index])
    # # pass
    # normal_vector = norm(points[9, :] - points[6, :])
    # z_ = -np.dot(normal_vector, vors.vertices[3])
    # plane = Plane(normal_vector=normal_vector, z_=z_)
    # fig = Axes3D(plt.figure())
    # plane.plot(fig)
    # index1 = [3, 2]
    # index2 = [9, 6]
    # fig.scatter(vors.vertices[index1, 0], vors.vertices[index1, 1], vors.vertices[index1, 2], marker='x')
    # fig.scatter(points[index2, 0], points[index2, 1], points[index2, 2], marker='o')
    # plt.show()
    bound = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    bound = Polyhedron(bound)
    points = np.random.rand(20, 3)
    fig = plt.figure()
    for i in range(6):
        all_points = points
        for plane in bound.planes:
            all_points = np.r_[all_points, mirror_projection(points, plane)]
        vor = Voronoi(all_points)
        polys = list()
        points = np.array([])
        for region in vor.regions:
            if region and -1 not in region:
                curr_points = vor.vertices[region]
                if np.min(curr_points) > -1e-3 and np.max(curr_points) < 1 + 1e-3:
                    polys.append(Polyhedron(curr_points))
                    points = np.r_[points, np.sum(polys[-1].points, 0) / polys[-1].points.shape[0]]
        points = points.reshape((-1, 3))
        ax = fig.add_subplot(2, 3, i + 1, projection='3d')
        for poly in polys:
            poly.plot(ax)
    plt.show()

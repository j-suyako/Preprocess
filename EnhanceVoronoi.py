from Lib.GeometryLib import *


class EhVornoi(Voronoi):
    pass


if __name__ == '__main__':
    np.random.seed(5)
    points = np.random.rand(10, 3)
    polys = [0] * points.shape[0]
    vors = Voronoi(points)
    # for pair in vors.ridge_dict:
    #     points_in_plane = vors.ridge_dict[pair]
    #     target_point_index = points_in_plane[0] if points_in_plane[0] != -1 else points_in_plane[1]
    #     this, that = pair
    #     normal_vector = norm(points[this, :] - points[that, :])
    #     z_ = -np.dot(normal_vector, points[target_point_index])
    # pass
    normal_vector = norm(points[9, :] - points[6, :])
    z_ = -np.dot(normal_vector, vors.vertices[3])
    plane = Plane(normal_vector=normal_vector, z_=z_)
    fig = Axes3D(plt.figure())
    plane.plot(fig)
    index1 = [3, 2]
    index2 = [9, 6]
    fig.scatter(vors.vertices[index1, 0], vors.vertices[index1, 1], vors.vertices[index1, 2], marker='x')
    fig.scatter(points[index2, 0], points[index2, 1], points[index2, 2], marker='o')
    plt.show()

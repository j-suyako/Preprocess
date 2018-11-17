import numpy as np


def norm(x):
    if len(x.shape) == 1:
        return x / np.sqrt(np.dot(x, x))
    else:
        return np.array([norm(e) for e in x])


def angel(vec1, vec2):
    if (np.abs(vec1) < 1e-5).all() or (np.abs(vec1) < 1e-5).all():
        return 0
    return np.arccos(np.sum(vec1 * vec2) / np.sqrt(np.sum(vec1 * vec1) * np.sum(vec2 * vec2)))


def rotate_z(points, theta):
    """点绕z轴逆时针旋转（从z轴正方向看）theta角

    :param points:
    :param theta:
    :return:
    """
    rotate_matrix = np.asarray([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    return np.dot(points, rotate_matrix)


def rotate_x(points, theta):
    rotate_matrix = np.asarray([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
    return np.dot(points, rotate_matrix)

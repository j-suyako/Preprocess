import numpy as np


def norm(x):
    if len(x.shape) == 1:
        if np.sum(np.abs(x)) <= 1e-5:
            return np.array([0, 0, 0])
        # x = np.array([0 if abs(e) < 1e-5 else e for e in x])
        x = x / np.max(np.abs(x))
        return x / np.sqrt(np.dot(x, x))
    else:
        return np.array([norm(e) for e in x])


def angel(vec1, vec2):
    if (np.abs(vec1) < 1e-5).all() or (np.abs(vec1) < 1e-5).all():
        return 0
    vec1 = norm(vec1)
    vec2 = norm(vec2)
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


def rotate(center: np.ndarray, points: np.ndarray, theta: iter):
    """rotate points by center

    :param center:
    :param points:
    :param theta: array_like
                iter of (theta_xy, theta_xz, theta_yz), means points first rotate theta_xy in xy plane,
                then rotate theta_xz in xz plane, finally rotate theta_yz in yz plane, counterclockwise
                is the positive direction
    :return:
    """
    # if len(points.shape) == 1:
    vec = np.mat(points - center).T
    theta_xy, theta_xz, theta_yz = theta
    a = np.mat([[np.cos(theta_xy), -np.sin(theta_xy), 0],
                [np.sin(theta_xy), np.cos(theta_xy), 0],
                [0, 0, 1]])
    b = np.mat([[np.cos(theta_xz), 0, -np.sin(theta_xz)],
                [0, 1, 0],
                [np.sin(theta_xz), 0, np.cos(theta_xz)]])
    c = np.mat([[1, 0, 0],
                [0, np.cos(theta_yz), -np.sin(theta_yz)],
                [0, np.sin(theta_yz), np.cos(theta_yz)]])
    return np.around(np.array(c * b * a * vec).T + center, 5)
    # else:
    #     return np.array([rotate(center, point, theta) for point in points])


def ridge_renum(mapping, ridge):
    if len(ridge) == 0:
        raise ValueError()
    ridge_map = dict(zip(mapping, list(range(len(mapping)))))
    if isinstance(ridge[0], list):
        return [[ridge_map[e] for e in r] for r in ridge]
    else:
        return [ridge_map[r] for r in ridge]


def continuous_block(a):
    res = []
    start, end = 0, 1
    while end < len(a):
        if a[end] - a[start] == end - start:
            pass
        elif end - start == 1:
            start = end
        else:
            res.append(a[start:end])
            start = end
        end += 1
    res.append(a[start:end])
    return res


if __name__ == "__main__":
    center = np.array([0, 0, 0])
    points = np.array([1, 0, 0])
    points_after_rotate = rotate(center, points, [np.pi / 2, 0, np.pi / 2])
    print(points_after_rotate)

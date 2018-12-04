import numpy as np


def transform(points: np.ndarray):
    if len(points.shape) == 1:
        return "[{}]".format(", ".join([str(e) for e in points]))
    else:
        return "[{}]".format(",\n ".join([transform(point) for point in points]))


if __name__ == "__main__":
    print(transform(np.array([[0, 0, 1],
                              [1, 1, 0]])))

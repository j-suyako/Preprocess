from config import LOGGER


def logging(process: str):
    def _logging(func):
        def __logging(*args, **kwargs):
            res = func(*args, **kwargs)
            LOGGER.info("{} finished.".format(process))
            return res
        return __logging
    return _logging

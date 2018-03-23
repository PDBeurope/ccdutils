
def str_is_int(i):
    try:
        int(i)
        return True
    except ValueError:
        return False


def str_to_int(i):
    try:
        return int(i)
    except ValueError:
        return 0


def str_to_float(f):
    try:
        return float(f)
    except ValueError:
        return 0.0

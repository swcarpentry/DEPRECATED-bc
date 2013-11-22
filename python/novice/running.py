def running(values):
    result = [0]
    for v in values[1:]:
        result.append(result[-1] + v)
    return result

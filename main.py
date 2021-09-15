import numpy as np


def _rand(p, a, b):
    return (p * a) % b


def _mean(lst):
    result = np.float64(0.0)
    for item in lst:
        result += item/len(lst)
    return result


def _std(lst, m=None):
    if m is None:
        m = _mean(lst)
    rate = np.float64(1.0) / (len(lst) - 1)
    result = np.float64(0.0)
    for item in lst:
        result += (item**2 - m**2) * rate
    return np.sqrt(result)


def f1(lst):
    max_seq_len = 0
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if lst[i] == lst[j]:
                if j - i > max_seq_len:
                    max_seq_len = j - i
                break
        break
        #if i >= len(lst) - max_seq_len:
        #    break
    return max_seq_len


def f2(lst):
    min_seq_len = len(lst)
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if lst[i] == lst[j]:
                if j - i < min_seq_len:
                    min_seq_len = j - i
                break
        break
    return min_seq_len


if __name__ == '__main__':
    history = [13]
    for i in range(100000):
        history.append(_rand(history[i], 73009, 63949))
    h_mean = _mean(history)
    h_std = _std(history, h_mean)
    h_max_seq_len = f1(history)
    h_min_seq_len = f2(history)
    print(f"Expectation: {h_mean:.5f} (effective mean)")
    print(f"Standard deviation: {h_std:.5f} (corrected)")
    print(f"Maximum unique sequence length: {h_max_seq_len}")
    print(f"Minimum unique sequence length: {h_min_seq_len}")

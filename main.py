import numpy as np
import matplotlib.pyplot as plt


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


def _hist(x, n_bins=20):
    n_bins = int(n_bins)
    if type(x) != np.ndarray:
        x = np.array(x)
    x_min = np.amin(x)
    x_range = np.amax(x) - x_min
    bin_step = float(x_range) / n_bins
    bins = [bin_step * i for i in range(n_bins + 1)]
    if bin_step > 0:
        result = [0] * n_bins
        for item in x:
            result[int((item - x_min) / bin_step)-1] += 1
    else:
        raise ZeroDivisionError
    return bins, result


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
        # 1 702551 958821
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
    # 13, 73009, 63949
    # 13 500933 200569 410819
    param_r0 = int(input(f"R0 parameter: "))
    print(f"(R * a) % b")
    param_a = int(input(f"a parameter: "))
    param_b = int(input(f"b parameter: "))
    history = [param_r0]
    for i in range(100000):
        history.append(_rand(history[i], param_a, param_b))
    h_mean = _mean(history)
    h_std = _std(history, h_mean)
    h_max_seq_len = f1(history)
    h_min_seq_len = f2(history)
    bins, hist = _hist(history)
    print(f"Expectation: {h_mean:.5f} (effective mean)")
    print(f"Standard deviation: {h_std:.5f} (corrected)")
    print(f"Maximum unique sequence length: {h_max_seq_len}")
    print(f"Minimum unique sequence length: {h_min_seq_len}")
    plt.figure().add_subplot()
    plt.gca().hist(bins[:-1], bins, weights=hist)
    plt.show()

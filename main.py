import queue

import numpy as np
import matplotlib.pyplot as plt


def _rand(p, a, b):
    return (p * a) % b


class SequenceAnalyzer:
    @staticmethod
    def mean(lst):
        result = np.float64(0.0)
        for i in range(len(lst)):
            result = result * (i / (i+1)) + (lst[i] / (i + 1))
        return result

    @staticmethod
    def std(lst, m=None):
        if m is None:
            m = SequenceAnalyzer.mean(lst)
        rate = np.float64(1.0) / (len(lst) - 1)
        result = np.float64(0.0)
        for item in lst:
            result += (item**2 - m**2) * rate
        return np.sqrt(result)

    @staticmethod
    def hist(x, n_bins=20):
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

    @staticmethod
    def lehmer_analysis(lst):
        # 1 702551 958821
        # 1 125566 276128
        # 8628 8632

        period = len(lst)
        max_unique = -1
        i = len(lst) - 1
        while i > 0:
            j = i - 1
            low_bound = max(0, i - period)
            while j >= low_bound:
                if lst[i] == lst[j]:
                    if i - j < period:
                        period = i - j
                    break
                j -= 1
            else:
                i -= 1
                continue
            break

        if period < len(lst):
            i = 0
            while lst[i] != lst[i+period]:
                i += 1
            else:
                max_unique = i + period
                if i > 0:
                    max_unique -= 1
        return period, max_unique

    @staticmethod
    def general_analysis(lst):
        periodical_ranges = []

        max_seq_len = -1
        min_period = len(lst)
        i = len(lst) - 1
        while i > 0:
            j = i - 1
            low_bound = max(0, i - min_period)
            while j >= low_bound:
                if lst[i] == lst[j]:
                    if i - j < min_period:
                        min_period = i - j
                    k = i
                    while True:
                        i -= 1
                        j -= 1
                        if not (j > 0 and lst[i] == lst[j]):
                            break
                    periodical_ranges.append(tuple((k, i-j, k-i)))
                    if max_seq_len < k - i:
                        max_seq_len = k - i
                    break
                j -= 1
            else:
                i -= 1

        def range_contains(x, p):
            return x[0] <= p < (x[0] + x[1])

        # ... range[0] + range[1]
        # range[0] + range[2] ...
        # ... range0[0] ... range1[0] + range1[1] if range0[2] < range0[1]
        sorted_ranges = []
        for item in periodical_ranges:
            sorted_ranges.append(tuple((item[0] - item[1] - item[2], item[1], item[2])))
        sorted_ranges = sorted(sorted_ranges, key=lambda item: item[0])

        unique_seq_lens = []
        start_points = queue.Queue(len(sorted_ranges))
        start_points.put(0)
        start_range = 0
        try:
            while True:
                point = start_points.get(False)
                while point > sorted_ranges[start_range][0] and \
                        not range_contains(sorted_ranges[start_range], point):
                    start_range += 1
                    if start_range >= len(sorted_ranges):
                        break
                else:
                    current = sorted_ranges[start_range]
                    if sorted_ranges[start_range][1] > sorted_ranges[start_range][2]:
                        high_bound = current[0] + current[1]
                        for i in range(start_range + 1, len(sorted_ranges)):
                            if range_contains(current, sorted_ranges[i][0]):
                                if range_contains(current, sorted_ranges[i][0] + sorted_ranges[i][1]):
                                    if sorted_ranges[i][0] + sorted_ranges[i][1] < high_bound:
                                        high_bound = sorted_ranges[i][0] + sorted_ranges[i][1]
                                    start_points.put(sorted_ranges[i][0] + sorted_ranges[i][2])
                            else:
                                break
                        unique_seq_lens.append(high_bound - point)
                        if high_bound == current[0] + current[1]:
                            start_points.put(high_bound)
                    else:
                        unique_seq_lens.append(current[0] + current[1] - point)
                        start_points.put(current[0] + current[2])
                    continue
                unique_seq_lens.append(len(lst) - point)
                break
        except queue.Empty:
            pass

        max_unique = max(unique_seq_lens)
        return min_period, max_seq_len, max_unique

    @staticmethod
    def geometry_test(x):
        counter = 0
        for i in range(0, len(x)-1, 2):
            if (x[i]**2 + x[i+1]**2) < 1:
                counter += 1
        return (2 * counter) / len(x)


def lab():
    # params = [13, 500933, 410819]
    params = [1, 125566, 276128]
    _labels = [f"R0 parameter: ", f"k parameter: ", f"b parameter: "]
    print(f"(R * k) % b")
    for i in range(len(params)):
        try:
            temp = int(input(_labels[i]))
            params[i] = temp if temp != 0 else params[i]
        except ValueError:
            continue

    seq_len = min(int(1e6), params[2])
    history = [params[0]]
    for i in range(seq_len):
        history.append(_rand(history[i], params[1], params[2]))

    # h_min_period, h_max_seq_len, h_max_unique_len = SequenceAnalyzer.general_analysis(history)
    # print(f"Maximum periodical sequence length: {h_max_seq_len}")
    # print(f"Minimum period: {h_min_period}")
    # print(f"Maximum not periodical sequence length: {h_max_unique_len}")

    h_period, h_unique_len = SequenceAnalyzer.lehmer_analysis(history)
    print(f"Period: {h_period}")
    print(f"Maximum not periodical sequence length: {h_unique_len}")

    for i in range(len(history)):
        history[i] /= params[2]
    h_mean = SequenceAnalyzer.mean(history)
    h_std = SequenceAnalyzer.std(history, h_mean)
    h_circle_hit = SequenceAnalyzer.geometry_test(history)
    bins, hist = SequenceAnalyzer.hist(history)
    print(f"Expectation: {h_mean:.5f} (effective mean)")
    print(f"Standard deviation: {h_std:.5f} (corrected)")
    print(f"Circle hit rate: {h_circle_hit:.5f}, pi/4={np.pi / 4}")
    plt.figure().add_subplot()
    plt.gca().hist(bins[:-1], bins, weights=hist)
    plt.show()


if __name__ == '__main__':
    lab()

"""while lst[i] == lst[j] and j > 0:
    i -= 1
    j -= 1"""
"""max_seq_len = 0
        max_period = 0
        i = len(lst) - 1
        while i > max_period:
            j = i - max_period - 1
            while j >= 0:
                if lst[i] == lst[j]:
                    if i - j > max_period:
                        max_period = i - j
                        k = max_seq_len
                        while k < max_period:
                            if lst[i] != lst[j]:
                                break
                            i -= 1
                            j -= 1
                        else:
                            i = j
                        max_seq_len = k
                        break
                j -= 1
            else:
                i -= 1"""

"""@staticmethod
def f2(lst):
    max_seq_len = 0
    min_period = len(lst)
    max_unique = -1
    i = len(lst) - 1
    sequence_start = i
    while i > 0:
        j = i - 1
        low_bound = max(0, i - min_period)
        while j >= low_bound:
            if lst[i] == lst[j]:
                if sequence_start - j < max_unique:
                    if sequence_start - j < max_unique - (sequence_start - i):
                        sequence_start = i
                    max_unique = max(sequence_start - j, max_unique - (sequence_start - i))
                if i - j > max_unique:
                    max_unique = i - j
                if i - j < min_period:
                    min_period = i - j
                k = i
                while True:
                    i -= 1
                    j -= 1
                    if not (j > 0 and lst[i] == lst[j]):
                        break
                if max_seq_len < k - i:
                    max_seq_len = k - i
                break
            j -= 1
        else:
            i -= 1
    return min_period, max_seq_len, max_unique"""

#!/usr/bin/env python3
"""Manhattan Tourist Problem (Python)"""


def manhattan_tourist(n, m, down, right):
    s = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(
                s[i - 1][j] + down[i - 1][j],
                s[i][j - 1] + right[i][j - 1]
            )

    return s[n][m]


if __name__ == "__main__":
    down = [
        [3, 2, 4, 0, 1],
        [3, 2, 4, 2, 0],
        [0, 7, 3, 4, 4],
        [3, 3, 0, 2, 1],
    ]
    right = [
        [1, 4, 5, 2],
        [0, 6, 5, 8],
        [4, 2, 5, 3],
        [1, 3, 2, 1],
        [3, 1, 1, 3],
    ]
    result = manhattan_tourist(4, 4, down, right)
    print(f"Maximum score path: {result}")
    print("This is the foundation of sequence alignment algorithms.")
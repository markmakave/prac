import sys
import numpy as np

np.random.seed(42)

m = int(sys.argv[1])
n = int(sys.argv[2])

A = np.random.randint(-100, 100, size=(m * n + 2), dtype=np.int32)
A[0] = m
A[1] = n

B = np.random.randint(-100, 100, size=(n + 1), dtype=np.int32)
B[0] = n

A.tofile('A.bin')
B.tofile('B.bin')

A = A[2:].reshape(m, n)
B = B[1:]

REF = np.dot(A, B)
REF = np.insert(REF, 0, REF.shape[0])
REF.tofile('REF.bin')



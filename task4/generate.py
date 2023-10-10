import sys
import numpy as np

np.random.seed(42)

n = int(sys.argv[1])

A = np.random.randint(-100, 100, size=(n*n + 1), dtype=np.int32)
A[0] = n
B = np.random.randint(-100, 100 + 1, size=(n*n + 1), dtype=np.int32)
B[0] = n

A.tofile('A.bin')
B.tofile('B.bin')

A = A[1:].reshape(n, n)
B = B[1:].reshape(n, n)

C = np.dot(A, B)
C = C.reshape(n*n)
C = np.insert(C, 0, n)
C.tofile('REF.bin')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import time
from PIL import Image

start = time.perf_counter()

A = 4
N = 21256             # Grid size = (N * A)^2
cropRGB = np.zeros((N, N, 3), dtype = np.uint8)
graphRGB = np.memmap('/kaggle/temp/graphRGB.dat', dtype = 'uint8', mode = 'w+', shape = (A * N, A * N, 3))
graphRGB.flush()

libFractalGPU = np.ctypeslib.load_library("libFractalGPU.so", ".")

fractalRootCalcGPU = libFractalGPU.fractalRootCalc
fractalRootCalcGPU.argtypes = [np.ctypeslib.ndpointer(), np.ctypeslib.ctypes.c_int, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double]
fractalOrderCalcGPU = libFractalGPU.fractalOrderCalc
fractalOrderCalcGPU.argtypes = [np.ctypeslib.ndpointer(), np.ctypeslib.ctypes.c_int, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double, np.ctypeslib.ctypes.c_double]

for i in range(A) :
    for j in range(A) :
        fractalRootCalcGPU(cropRGB, N, -2.0 + i * 4.0 / A, -2.0 + (i + 1) * 4.0 / A, -2.0 + j * 4.0 / A, -2.0 + (j + 1) * 4.0 / A)
        graphRGB[(j * N) : ((j + 1) * N), (i * N) : ((i + 1) * N)] = cropRGB
        graphRGB.flush()

Image.fromarray(graphRGB).save('rootRGB.png')
print("rootRGB.png written.")

for i in range(A) :
    for j in range(A) :
        fractalOrderCalcGPU(cropRGB, N, -2.0 + i * 4.0 / A, -2.0 + (i + 1) * 4.0 / A, -2.0 + j * 4.0 / A, -2.0 + (j + 1) * 4.0 / A)
        graphRGB[(j * N) : ((j + 1) * N), (i * N) : ((i + 1) * N)] = cropRGB
        graphRGB.flush()

Image.fromarray(graphRGB).save('orderRGB.png')
print("orderRGB.png written.")

end = time.perf_counter()
print(f"Time: {end - start:.5f}s")

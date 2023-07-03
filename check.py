import ctypes

lib = ctypes.CDLL('./libexample.so')
result = lib.grating_coupler_period(float(1550), float(1.44), float(1.5), float(60), float(1))
print(result)
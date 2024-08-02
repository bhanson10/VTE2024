import ctypes as ct

lib = ct.CDLL("./hello_world.so")
lib.hello_world()
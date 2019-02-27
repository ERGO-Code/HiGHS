import ctypes
import os
highs = ctypes.cdll.LoadLibrary("libhighs.so") # highs lib folder must be in "LD_LIBRARY_PATH" environment variable

highs.callhighs.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double))

def call_highs(colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue):
   global highs
   n_col = len(colcost)
   n_row = len(rowlower)
   n_nz = len(aindex)

   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_astart = ctypes.c_int * (n_col + 1)
   int_array_type_aindex = ctypes.c_int * n_nz
   dbl_array_type_avalue = ctypes.c_double * n_nz

   highs.callhighs(ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), dbl_array_type_col(*colcost), dbl_array_type_col(*collower), dbl_array_type_col(*colupper), dbl_array_type_row(*rowlower), dbl_array_type_row(*rowupper), int_array_type_astart(*astart), int_array_type_aindex(*aindex), dbl_array_type_avalue(*avalue))

import ctypes
import os

highslib = ctypes.cdll.LoadLibrary("libhighs.so") # highs lib folder must be in "LD_LIBRARY_PATH" environment variable

highslib.Highs_mipCall.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int))
highslib.Highs_call.restype = ctypes.c_int

def Highs_mipCall(colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, integrality):
   global highslib
   n_col = len(colcost)
   n_row = len(rowlower)
   n_nz = len(aindex)
   rowwise = 0
   
   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_astart = ctypes.c_int * n_col
   int_array_type_aindex = ctypes.c_int * n_nz
   dbl_array_type_avalue = ctypes.c_double * n_nz

   int_array_type_col = ctypes.c_int * n_col
   int_array_type_row = ctypes.c_int * n_row

   col_value = [0] * n_col
   col_dual = [0] * n_col

   row_value = [0] * n_row
   row_dual = [0] * n_row

   col_basis = [0] * n_col
   row_basis = [0] * n_row

   return_val = 0

   col_value = dbl_array_type_col(*col_value)
   row_value = dbl_array_type_row(*row_value)

   retcode = highslib.Highs_mipCall(
       ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), ctypes.c_int(rowwise), 
       dbl_array_type_col(*colcost), dbl_array_type_col(*collower), dbl_array_type_col(*colupper), 
       dbl_array_type_row(*rowlower), dbl_array_type_row(*rowupper), 
       int_array_type_astart(*astart), int_array_type_aindex(*aindex), dbl_array_type_avalue(*avalue),
       int_array_type_col(*integrality), 
       col_value, row_value, ctypes.byref(ctypes.c_int(return_val)))
   return retcode, list(col_value), list(row_value)


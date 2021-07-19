import ctypes
import os
from ctypes.util import find_library

highslib = ctypes.cdll.LoadLibrary(ctypes.util.find_library("highs"))

# highs lib folder must be in "LD_LIBRARY_PATH" environment variable
# ============
# Highs_lpCall
highslib.Highs_lpCall.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, 
                                  ctypes.c_int, ctypes.c_double,
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
highslib.Highs_lpCall.restype = ctypes.c_int

def Highs_lpCall(col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value):
   global highslib
   n_col = len(col_cost)
   n_row = len(row_lower)
   n_nz = len(a_index)
   a_format = 1
   sense = 1
   offset = 0
   
   # In case a_start has the fictitious start of column n_col
   a_start_length = len(a_start)

   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_a_start = ctypes.c_int * a_start_length
   int_array_type_a_index = ctypes.c_int * n_nz
   dbl_array_type_a_value = ctypes.c_double * n_nz

   int_array_type_col = ctypes.c_int * n_col
   int_array_type_row = ctypes.c_int * n_row

   col_value = [0] * n_col
   col_dual = [0] * n_col

   row_value = [0] * n_row
   row_dual = [0] * n_row

   col_basis = [0] * n_col
   row_basis = [0] * n_row

   model_status = ctypes.c_int(0)

   col_value = dbl_array_type_col(*col_value)
   col_dual = dbl_array_type_col(*col_dual)
   row_value = dbl_array_type_row(*row_value)
   row_dual = dbl_array_type_row(*row_dual)
   col_basis = int_array_type_col(*col_basis)
   row_basis = int_array_type_row(*row_basis)

   return_status = highslib.Highs_lpCall(
      ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), ctypes.c_int(a_format),
      ctypes.c_int(sense), ctypes.c_double(offset),
      dbl_array_type_col(*col_cost), dbl_array_type_col(*col_lower), dbl_array_type_col(*col_upper), 
      dbl_array_type_row(*row_lower), dbl_array_type_row(*row_upper), 
      int_array_type_a_start(*a_start), int_array_type_a_index(*a_index), dbl_array_type_a_value(*a_value),
      col_value, col_dual, 
      row_value, row_dual, 
      col_basis, row_basis, ctypes.byref(model_status))
   return return_status, model_status.value, list(col_value), list(col_dual), list(row_value), list(row_dual), list(col_basis), list(row_basis)

# =============
# Highs_mipCall
highslib.Highs_mipCall.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, 
                                  ctypes.c_int, ctypes.c_double,
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int))
highslib.Highs_mipCall.restype = ctypes.c_int

def Highs_mipCall(col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality):
   global highslib
   n_col = len(col_cost)
   n_row = len(row_lower)
   n_nz = len(a_index)
   a_format = 1
   sense = 1
   offset = 0

   # In case a_start has the fictitious start of column n_col
   a_start_length = len(a_start)

   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_a_start = ctypes.c_int * a_start_length
   int_array_type_a_index = ctypes.c_int * n_nz
   dbl_array_type_a_value = ctypes.c_double * n_nz

   int_array_type_col = ctypes.c_int * n_col
   int_array_type_row = ctypes.c_int * n_row

   col_value = [0] * n_col
   col_dual = [0] * n_col

   row_value = [0] * n_row
   row_dual = [0] * n_row

   col_basis = [0] * n_col
   row_basis = [0] * n_row

   model_status = ctypes.c_int(0)

   col_value = dbl_array_type_col(*col_value)
   row_value = dbl_array_type_row(*row_value)

   return_status = highslib.Highs_mipCall(
      ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), ctypes.c_int(a_format), 
      ctypes.c_int(sense), ctypes.c_double(offset),
      dbl_array_type_col(*col_cost), dbl_array_type_col(*col_lower), dbl_array_type_col(*col_upper), 
      dbl_array_type_row(*row_lower), dbl_array_type_row(*row_upper), 
      int_array_type_a_start(*a_start), int_array_type_a_index(*a_index), dbl_array_type_a_value(*a_value),
      int_array_type_col(*integrality), 
      col_value, row_value, ctypes.byref(model_status))
   return return_status, model_status.value, list(col_value), list(row_value)

# ============
# Highs_qpCall
highslib.Highs_qpCall.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                  ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
highslib.Highs_call.restype = ctypes.c_int

def Highs_qpCall(col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value):
   global highslib
   n_col = len(col_cost)
   n_row = len(row_lower)
   n_nz = len(a_index)
   q_n_nz = len(q_index)
   a_format = 1
   q_format = 1
   sense = 1
   offset = 0
   
   # In case a_start or q_start has the fictitious start of column n_col
   a_start_length = len(a_start)
   q_start_length = len(q_start)

   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_a_start = ctypes.c_int * a_start_length
   int_array_type_a_index = ctypes.c_int * n_nz
   dbl_array_type_a_value = ctypes.c_double * n_nz

   int_array_type_q_start = ctypes.c_int * q_start_length
   int_array_type_q_index = ctypes.c_int * q_n_nz
   dbl_array_type_q_value = ctypes.c_double * q_n_nz

   int_array_type_col = ctypes.c_int * n_col
   int_array_type_row = ctypes.c_int * n_row

   col_value = [0] * n_col
   col_dual = [0] * n_col

   row_value = [0] * n_row
   row_dual = [0] * n_row

   col_basis = [0] * n_col
   row_basis = [0] * n_row

   model_status = ctypes.c_int(0)

   col_value = dbl_array_type_col(*col_value)
   col_dual = dbl_array_type_col(*col_dual)
   row_value = dbl_array_type_row(*row_value)
   row_dual = dbl_array_type_row(*row_dual)
   col_basis = int_array_type_col(*col_basis)
   row_basis = int_array_type_row(*row_basis)

   return_status = highslib.Highs_qpCall(
      ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), ctypes.c_int(q_n_nz), 
      ctypes.c_int(a_format), ctypes.c_int(q_format),
      ctypes.c_int(sense), ctypes.c_double(offset),
      dbl_array_type_col(*col_cost), dbl_array_type_col(*col_lower), dbl_array_type_col(*col_upper), 
      dbl_array_type_row(*row_lower), dbl_array_type_row(*row_upper), 
      int_array_type_a_start(*a_start), int_array_type_a_index(*a_index), dbl_array_type_a_value(*a_value),
      int_array_type_q_start(*q_start), int_array_type_q_index(*q_index), dbl_array_type_q_value(*q_value),
      col_value, col_dual, 
      row_value, row_dual, 
      col_basis, row_basis, ctypes.byref(model_status))
   return return_status, model_status.value, list(col_value), list(col_dual), list(row_value), list(row_dual), list(col_basis), list(row_basis)

# ==========================================
# Highs_lpDimMpsRead and Highs_lpDataMpsRead
highslib.Highs_lpDimMpsRead.argtypes = (ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int),
                                        ctypes.POINTER(ctypes.c_int))
highslib.Highs_lpDimMpsRead.restype = ctypes.c_int

highslib.Highs_lpDataMpsRead.argtypes = (ctypes.c_int, ctypes.c_int, 
                                         ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                         ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_double))
highslib.Highs_lpDataMpsRead.restype = ctypes.c_int

def Highs_lpDimMpsRead():
   global highslib
   n_col = ctypes.c_int(0)
   n_row = ctypes.c_int(0)
   n_nz = ctypes.c_int(0)

   return_status = highslib.Highs_lpDimMpsRead(ctypes.byref(n_col),
                                               ctypes.byref(n_row),
                                               ctypes.byref(n_nz))
   return return_status, n_col.value, n_row.value, n_nz.value

def Highs_lpMpsRead():
   global highslib
   status, n_col, n_row, n_nz = Highs_lpDimMpsRead() 

   sense = 1
   offset = 0
   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_a_start = ctypes.c_int * (n_col+1)
   int_array_type_a_index = ctypes.c_int * n_nz
   dbl_array_type_a_value = ctypes.c_double * n_nz
   col_cost = [0] * n_col
   col_lower = [0] * n_col
   col_upper = [0] * n_col
   row_lower = [0] * n_row
   row_upper = [0] * n_row
   a_start = [0] * (n_col+1)
   a_index = [0] * n_nz
   a_value = [0] * n_nz
   col_cost = dbl_array_type_col(*col_cost)
   col_lower = dbl_array_type_col(*col_lower)
   col_upper = dbl_array_type_col(*col_upper)
   row_lower = dbl_array_type_row(*row_lower)
   row_upper = dbl_array_type_row(*row_upper)
   a_start = int_array_type_a_start(*a_start)
   a_index = int_array_type_a_index(*a_index)
   a_value = dbl_array_type_a_value(*a_value)

   return_status = highslib.Highs_lpDataMpsRead(ctypes.c_int(n_col), ctypes.c_int(n_row),
                                          ctypes.byref(ctypes.c_int(sense)), ctypes.byref(ctypes.c_double(offset)),
                                          col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value);
   return return_status, n_col, n_row, sense, offset, list(col_cost), list(col_lower), list(col_upper), list(row_lower), list(row_upper), list(a_start), list(a_index), list(a_value)


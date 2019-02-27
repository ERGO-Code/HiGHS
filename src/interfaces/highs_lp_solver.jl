
function call_highs(cc, cl, cu, rl, ru, astart, aindex, avalue)
   n_col = convert(Int32, size(cc, 1))
   n_row = convert(Int32, size(rl, 1))
   n_nz = convert(Int32, size(aindex, 1))

   colcost = convert(Array{Cdouble}, cc)
   collower = convert(Array{Cdouble}, cl)
   colupper = convert(Array{Cdouble}, cu)

   rowlower = convert(Array{Cdouble}, rl)
   rowupper = convert(Array{Cdouble}, ru)
   matstart = convert(Array{Int32}, astart)
   matindex = convert(Array{Int32}, aindex)
   matvalue = convert(Array{Cdouble}, avalue)

   ccall((:callhighs, "libhighs.so"), Cvoid, (Int32, Int32, Int32, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Int32},Ptr{Int32},Ptr{Cdouble}),
   n_col, n_row, n_nz, colcost, collower, colupper, rowlower, rowupper, matstart, matindex, matvalue)
end
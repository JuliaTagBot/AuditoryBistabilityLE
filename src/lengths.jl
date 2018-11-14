export percept_lengths

function findlengths(x)
  indices = [1; findall(diff(x) .!= 0) .+ 1; length(x)+1]
  filter(!iszero,diff(indices)), x[indices[1:end-1]]
end

mergelengths(threshold) = lens -> mergelengths(lens,threshold)
function mergelengths(lens,vals,threshold)
  mergedlen = similar(lens)
  mergedval = similar(vals)

  i = 0
  for (len,val) in zip(lens,vals)
    if i > 0 && (len < threshold || mergedval[i] == val)
      mergedlen[i] += len
    else
      i += 1
      mergedlen[i] = len
      mergedval[i] = val
    end
  end

  mergedlen[1:i], mergedval[1:i]
end

# interpolate_times(x::MetaUnion{AxisArray};to) = 
#   interpolate_times(x,times(x),times(to),dim=axisdim(x,Axis{:time}))

# function interpolate_times(x::AbstractArray, times::AbstractVector,
#                            to::AbstractVector; dim=ndims(x))
#   @assert size(x)[dim] == length(times)
#   y = similar(x,float(eltype(x)),(size(x)[1:dim-1]...,length(to),
#                                   size(x)[dim:end]...))
#   x1 = 1
#   for i in 1:size(y,dim)
#     while x1 <= length(times) && times[x1] < to[i]; x1+=1; end
#     if x1 == 1
#       selectdim(y,dim,i) .= selectdim(x,dim,1)
#     elseif x1 > length(times)
#       selectdim(y,dim,i) .= selectdim(x,dim,size(x,dim))
#     else
#       a,b = selectdim(x,dim,x1-1), selectdim(x,dim,x1)
#       k,h = to[i] - times[x1-1], times[x1] - to[i]
#       selectdim(y,dim,i) .= a.*(h/(k+h)) .+ b.*(k/(k+h))
#     end
#   end
#   y
# end

function percept_lengths(counts; threshold = 0.45, 
                         min_length_ms=750, min_length=min_length_ms*ms)
  lens,vals = findlengths(Array(counts .< threshold))
  slens = lens * ustrip(uconvert(s,Δt(counts)))

  mergelengths(slens,vals,ustrip(uconvert(s,min_length)))
end


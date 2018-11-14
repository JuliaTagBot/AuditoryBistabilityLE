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

function percept_lengths(counts; threshold = 0.45, 
                         min_length_ms=750, min_length=min_length_ms*ms)
  lens,vals = findlengths(Array(counts .< threshold))
  slens = lens * ustrip(uconvert(s,Δt(counts)))

  mergelengths(slens,vals,ustrip(uconvert(s,min_length)))
end


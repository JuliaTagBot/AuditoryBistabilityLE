using CodecZlib
using Statistics
export compress, compress!, decompress

struct CompressedMask{T}
  size::Tuple{Int,Int}
  times::T
  data::Array{UInt8}
end

compress(x::AbstractMatrix) = compress!(deepcopy(x))

function compress!(x::AbstractMatrix)
  x ./= maximum(x)
  sortedx = sort(vec(x))
  x[x .< 0] .= 0
  quantized = floor.(UInt8,x.*typemax(UInt8))
  CompressedMask(size(quantized), times(x),
                 transcode(ZlibCompressor,vec(quantized)))
end

function decompress(x::CompressedMask)
  quantized = transcode(ZlibDecompressor,x.data)
  mask = quantized ./ typemax(UInt8)
  AxisArray(reshape(mask,x.size...), Axis{:time}(x.times))
end



  

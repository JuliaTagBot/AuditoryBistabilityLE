using CodecZlib
using Statistics
export compress, compress!, decompress

struct CompressedMask{A}
  size::Tuple{Int,Int}
  axes::A
  data::Array{UInt8}
end

compress_axes(x::MetaUnion{AxisArray}) = AxisArrays.axes(x)
compress_axes(x) = nothing

function compress(x::AbstractMatrix)
  quantized = Array{UInt8}(undef,prod(size(x)))
  mx = maximum(x)
  for (i,ii) in enumerate(eachindex(x))
    quantized[i] = floor(UInt,max(0,x[ii]/mx)*typemax(UInt8))
  end
  CompressedMask(size(x), compress_axes(x),
                 transcode(ZlibCompressor,quantized))
end

withaxes(x,axes) = AxisArray(x, axes...)
withaxes(x,::Nothing) = x

function decompress(x::CompressedMask)
  mask = transcode(ZlibDecompressor,x.data) ./ typemax(UInt8)
  withaxes(reshape(mask,x.size...), x.axes)
end
  

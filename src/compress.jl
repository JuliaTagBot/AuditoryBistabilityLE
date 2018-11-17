using CodecZlib
using Statistics
export compress, decompress

struct CompressedMask{A}
  size::Tuple{Int,Int}
  axes::A
  data::Array{UInt8}
end

compress_axes(x::MetaUnion{AxisArray}) = AxisArrays.axes(x)
compress_axes(x) = nothing

function compress(x::AbstractMatrix)
  quantized = Array{UInt8}(undef,size(x))
  @. quantized = floor(UInt8,max(0,x/$maximum(x))*typemax(UInt8))
  CompressedMask(size(x), compress_axes(x),
                 transcode(ZlibCompressor,vec(quantized)))
end

withaxes(x,axes) = AxisArray(x, axes...)
withaxes(x,::Nothing) = x

function decompress(x::CompressedMask)
  mask = transcode(ZlibDecompressor,x.data) ./ typemax(UInt8)
  withaxes(reshape(mask,x.size...), x.axes)
end
  
function ShammaModel.audiospect(x::CompressedMask,settings)
  settings = read_settings(settings)
  audiospect(decompress(x);settings.freqs.analyze...)
end


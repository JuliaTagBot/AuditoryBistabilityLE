using CodecZlib
using Statistics
export compress, compress!, uncompress

struct CompressedMask
  size::Tuple{Int,Int}
  data::Array{UInt8}
end

compress(x::AbstractMatrix) = compress!(deepcopy(x))

function compress!(x::AbstractMatrix)
  x ./= maximum(x)
  sortedx = sort(vec(x))
  minthresh = quantile(sortedx,0.5,sorted=true)
  maxthresh = quantile(sortedx,0.85,sorted=true)
  x[x .< clamp(0.1,minthresh,maxthresh)] .= 0
  quantized = floor.(UInt8,x.*typemax(UInt8))
  CompressedMask(size(quantized), transcode(ZlibCompressor,vec(quantized)))
end

function uncompress(x::CompressedMask)
  quantized = transcode(ZlibDecompressor,x.data)
  mask = quantized ./ typemax(UInt8)
  reshape(mask,x.size...)
end



  

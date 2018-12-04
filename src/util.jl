read_settings(x) = as_namedtuple(x)
read_settings(x::String) = as_namedtuple(TOML.parsefile(x))
as_namedtuple(xs) = xs
function as_namedtuple(xs::Dict{<:AbstractString,<:Any})
  kt = ((Symbol(x) for x in keys(xs))...,)
  vt = (values(xs)...,)
  NamedTuple{(kt...,)}(as_namedtuple.(vt))
end

read_params(x) = x
function read_params(x::DataFrame)
  @assert size(x,1) == 1 "Only a single row of parameters can be provided."
  Dict(k => x[k][1] for k in names(x))
end


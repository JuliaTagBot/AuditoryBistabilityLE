using ShammaModel
using AxisArrays
const reasonable_response_maximum = 100
export percept_lengths, bandwidth_ratio

function bound(x,min,max)
  y = 1/(1+exp(-8((x - min)/(max - min) - 0.5)))
  miny = 1/(1+exp(-8(-0.5)))
  maxy = 1/(1+exp(-8(0.5)))
  min + (max-min)*(clamp(y,miny,maxy) - miny)/(maxy - miny)
end

function percept_lengths(result::NamedTuple,settings)
  settings = read_settings(settings)
  ratio = bandwidth_ratio(result,settings)
  percept_lengths(ratio,settings)
end

percept_lengths(ratio::AbstractVector,settings) =
  percept_lengths(ratio;settings.percept_lengths...)

function fixlength(x,n;fill=0)
  if length(x) > n
    view(x,1:n)
  elseif length(x) < n
    append!(convert(Array{Union{eltype(x),Missing}},x),
            Base.fill(fill,n-length(x)))
  else
    x
  end
end

function dominant_channels(sp;level_threshold=0.9,threshold=0.25,
                           window=500ms,delta=250ms,thresh_ratio=1.0,
                            max_channels=10)
  x = map_windowing(sp,length=window,step=delta,flatten=true) do window
    thresh_level = thresh_ratio*quantile(vec(sp),1 - threshold/size(sp,2))
    # @show thresh_level
    levels = dropdims(mapslices(x -> quantile(x,level_threshold),
                                Array(window),dims=1), dims=1)
    # @show sort(levels,rev=true)[1:min(end,20)]
    # thresh_level = thresh_ratio*
    #   quantile(levels,1 - threshold/length(levels))

    x = fixlength(findall(levels .> thresh_level),max_channels,fill=missing)
    # @show x
    x
  end 
end

function bandwidth_ratio(spmask::AbstractMatrix, sp::AbstractMatrix,
                         settings)
  settings = read_settings(settings)
  bandwidth_ratio(spmask, sp; settings.bandwidth_ratio...)
end

function percept_lengths(spmask::AbstractMatrix, sp::AbstractMatrix,
                         settings)
  ratio = bandwidth_ratio(spmask, sp, settings)
  percept_lengths(ratio, settings)
end

function bandwidth_ratio(spmask, sp; threshold=1.5,
                         window_ms=500,window=window_ms*ms,
                         full_band_ratio=2,
                         level_threshold=0.9,
                         thresh_ratio=1.0,
                         delta_ms=250,delta=delta_ms*ms)
  @assert freqs(spmask) == freqs(sp) "Frequency axes are not equal."
  @assert times(spmask) == times(sp) "Time axes are not equal."

  println("fullchan")
  fullchan = dominant_channels(sp,window=window*full_band_ratio,
                               delta=delta,threshold=threshold,
                               level_threshold=level_threshold,
                               thresh_ratio=thresh_ratio)

  println("maskchan")
  maskchan = dominant_channels(spmask,window=window,delta=delta,
                               threshold=threshold,
                               level_threshold=level_threshold,
                               thresh_ratio=thresh_ratio)

  dim = axisdim(fullchan,Axis{:time})
  others = setdiff(1:ndims(fullchan),[dim])
  fullband = mapslices(Array(fullchan),dims=others) do x
    maximum(skipmissing(x)) - minimum(skipmissing(x))
  end |> x -> AxisArray(vec(x),Axis{:time}(times(fullchan)))

  maskband = AxisArray(Array{Float64}(undef,size(maskchan,dim)),
                       Axis{:time}(times(maskchan)))
  ratio = AxisArray(Array{Float64}(undef,size(maskchan,dim)),
                    Axis{:time}(times(maskchan)))
  winwidth = window*full_band_ratio
  for (i,t) in enumerate(times(maskchan))
    val = (-window*full_band_ratio*0.51 .. window*full_band_ratio*0.51) + t
    maskwin = view(maskchan,Axis{:time}(t))
    fullwin = view(fullchan,Axis{:time}(val))
    maskcount = skipmissing(intersect(maskwin,fullwin))
    if isempty(maskcount)
      maskband[i] = 0
    else
      maskband[i] = maximum(maskcount) - minimum(maskcount)
    end
    ratio[i] = maskband[i] / maximum(view(fullband,Axis{:time}(val)))
  end

  ratio, fullband, maskband
end

function meanabs(A,n)
  result = reduce((x,y) -> x+abs(y),A,dims=n,init=zero(real(eltype(A))))
  result ./= size(A,n)
  dropdims(result,dims=n)
end

function findweights(condition,x,normalize_start)
  weights = if condition == :scales
    AxisArray(meanabs(x,axisdim(x,Axis{:freq})) .*
              ustrip.(uconvert.(cycoct,ShammaModel.scales(x)))',
              AxisArrays.axes(x,Axis{:time}), AxisArrays.axes(x,Axis{:scale}))
  elseif condition ∈ [:freqs,:track]
    x
  end

  if normalize_start > 0s
    after_buildup = normalize_start .. last(times(weights))
    weights .-= minimum(weights[after_buildup])
    weights ./= maximum(weights[after_buildup])
  end
  weights
end

function remove_key_prefix!(prefix,dict)
  for key in keys(dict)
    if startswith(string(key),prefix)
      dict[Symbol(string(key)[length(prefix)+1:end])] = dict[key]
    end
  end
end

apply_bistable(x,args...;kwds...) = apply_bistable!(deepcopy(x),args...;kwds...)
function apply_bistable!(x,condition,params;
                         input_bound=(0.005,0.1),
                         lowpass=1.5,lowpass_order=3,
                         interactive=false,
                         intermediate_results=interactive,
                         normalize_start_s=-1.0,
                         normalize_start=normalize_start_s*s,
                         progressbar=interactive)

  if condition == :freqs
    remove_key_prefix!("f_",params)
  elseif condition == :scales
    remove_key_prefix!("s_",params)
  elseif condition == :track
    remove_key_prefix!("t_",params)
  end

  if iszero(params[:c_a]) && iszero(params[:c_m]) && iszero(params[:c_σ])
    return (result=x,)
  end

  noise_params = Dict(
    :τ_σ => params[:τ_σ], :c_σ => params[:c_σ]
  )
  adapt_params = Dict(
    :c_m => params[:c_m], :τ_m => params[:τ_m],
    :c_a => params[:c_a], :τ_a => params[:τ_a],
    :c_x => params[:c_x], :τ_x => params[:τ_x]
  )

  weights = findweights(condition,x,normalize_start)
  weights .= bound.(weights,input_bound...)
  input_weights = copy(weights)

  wn = drift(weights;noise_params...,progressbar=progressbar)
  wna,a,m = adaptmi(wn;W_m=weighting(x,condition,params),
                    shape_y = x -> max(0,x),progressbar=progressbar,
                    adapt_params...)

  low = digitalfilter(Lowpass(lowpass;fs=ustrip(uconvert(Hz,1/Δt(wna)))),
                      Butterworth(lowpass_order))
  wna_low = filt!(similar(wna),low,wna)

  # shouldn't the below instead be replaced with:
  x .*= wna_low./weights

  if intermediate_results
    (result=x,inweights=input_weights,outweights=wna_low,adapt=a,inhibit=m)
  else
    (result=x,)
  end

end

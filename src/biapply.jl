using ShammaModel
using AxisArrays
const reasonable_response_maximum = 100

function bound(x,min,max)
  y = 1/(1+exp(-8((x - min)/(max - min) - 0.5)))
  miny = 1/(1+exp(-8(-0.5)))
  maxy = 1/(1+exp(-8(0.5)))
  min + (max-min)*(clamp(y,miny,maxy) - miny)/(maxy - miny)
end

function strengths(tracks,tracks_lp;kwds...)
  strengths = map_components(tracks,tracks_lp;kwds...) do components
    component_means(components)
  end

  AxisArray(hcat(strengths...),AxisArrays.axes(strengths,Axis{:time}))
end

function ratio_to_lengths(ratio,threshold=2,min_length=1s)
  percept_lengths(AxisArray(ratios .< threshold,
                            AxisArrays.axes(ratios,Axis{:time})),min_length)
end

function component_ratio(tracks,tracks_lp;min_length=1s,
                         intermediate_results=false,kwds...)
  map_components(tracks,tracks_lp;kwds...) do components
    strengths = sort(component_means(components),rev=true)
    strengths[1] / sum(strengths[2:end])
  end
end

function estimate_bandwidth(sp;threshold=0.25,window=500ms,step=250ms)
  map_windowing(sp,length=window,step=step) do window
    means = dropdims(mean(window,dims=1),dims=1)
    peak = maximum(means)
    over = findall(means .> threshold*peak)
    maximum(over) - minimum(over) + 1
  end
end

function component_bandwidth_ratio(cs,tracks,tracks_lp;min_length=1s,
                                   threshold=0.25,
                                   progressbar=false,window=500ms,
                                   step=250ms)
  crmask = mask(cs,tracks,tracks_lp,window=window,step=step)
  spmask = audiospect(crmask,progressbar=progressbar)
  sp = audiospect(cs,progressbar=progressbar)

  fullband = estimate_bandwidth(sp,window=window,step=step,
                               threshold=threshold)
  maskband = estimate_bandwidth(spmask,window=window,step=step,
                               threshold=threshold)

  AxisArray(maskband ./ fullband,AxisArrays.axes(fullband,Axis{:time}))
end

function meanabs(A,n)
  result = reduce((x,y) -> x+abs(y),A,dims=n,init=zero(real(eltype(A))))
  result ./= size(A,n)
  dropdims(result,dims=n)
end

function findweights(condition,x)
  if condition == :scales
    AxisArray(meanabs(x,axisdim(x,Axis{:freq})) .*
              ustrip.(uconvert.(cycoct,ShammaModel.scales(x)))',
              AxisArrays.axes(x,Axis{:time}), AxisArrays.axes(x,Axis{:scale}))
  elseif condition ∈ [:freqs,:track]
    x
  end
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

  weights = findweights(condition,x)
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

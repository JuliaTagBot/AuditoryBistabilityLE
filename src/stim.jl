using Unitful
using SampledSignals

export aba_, ab, stimulus

const sr = 8000

silence(len) = SampleBuf(zeros(SampledSignals.inframes(Int,len,sr)),sr)
function tone(freq,len) 
  freq_ = SampledSignals.inHz(freq)
  len_ = SampledSignals.inseconds(len)
  SampleBuf(sin.(2π.*freq_.*range(0,stop=len_,step=1/sr)), sr)
end

function ramp(xs,len)
  n = SampledSignals.inframes(Int,len,samplerate(xs))
  t = range(0,stop=1,length=n)
  xs[1:n] .*= cos.(π.*(t.-1)./2)
  xs[end-n+1:end] .*= cos.(π.*t./2)
  xs
end

# TODO: add a test in sampled signals for .^2 and other power operators, as
# there seems to be a bug
normpower(x) = x ./ sqrt.(mean(x.*x,dims=1))
amplify(ratio) = x -> x.*ratio

function stimulus(total_len,repeats,freq,delta;tone_len_fraction=0.5,pattern="ab",
                  ramp_len_ms=0,ramp_len=ramp_len_ms*ms)
  @assert ramp_len isa Unitful.Time
  if pattern == "ab"
    ab(tone_len_fraction*total_len,(1-tone_len_fraction)*total_len,1,
       repeats,freq,delta;ramp_len=ramp_len)
  elseif pattern == "aba_"
    aba_(tone_len_fraction*total_len,(1-tone_len_fraction)*total_len,
         repeats,freq,delta;ramp_len=ramp_len)
  else
    error("Unexpected stimulus pattern '$pattern'")
  end
end

function aba_(tone_len,spacing_len,repeats,freq,delta;ramp_len=0ms)
  a_freq=freq
  b_freq=freq*(2^(delta/12))
  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  a = ramp_len > 0ms ? ramp(a,25ms) : a
  b = ramp_len > 0ms ? ramp(b,25ms) : b

  aba_ = [a;silence(spacing_len);
          b;silence(spacing_len);
          a;silence(2spacing_len+tone_len)]

  aba_seq = silence(0s)
  for i = 1:repeats
    aba_seq = [aba_seq; aba_]
  end

  aba_seq
end

function ab(tone_len,spacing_len,offset_ratio,repeats,freq,delta;ramp_len=0ms)
  @assert 0 <= offset_ratio <= 1

  a_freq=freq
  b_freq=freq*(2^(delta/12))
  a = tone(a_freq,tone_len)
  b = tone(b_freq,tone_len)
  a = ramp_len > 0ms ? ramp(a,25ms) : a
  b = ramp_len > 0ms ? ramp(b,25ms) : b

  offset = silence(offset_ratio*(tone_len+spacing_len))
  ab = [a;offset].+[offset;b]
  ab_len = 2(tone_len+spacing_len)
  ab = [ab; silence(ab_len - duration(ab)*s)]

  ab_seq = silence(0s)
  for i = 1:repeats
    ab_seq = [ab_seq; ab]
  end

  ab_seq
end

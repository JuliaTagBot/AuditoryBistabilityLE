using Unitful
using SampledSignals

export alter_ab, sync_ab, buildup_aba, aba, ab, ideal_ab

const sr = 8000

silence(len) = SampleBuf(zeros(SampledSignals.inframes(Int,len,sr)),sr)
function tone(freq,len) 
  freq_ = SampledSignals.inHz(freq)
  len_ = SampledSignals.inseconds(len)
  SampleBuf(sin.(2π.*freq_.*range(0,stop=len_,step=1/sr)), sr)
end

# TODO: add a test in sampled signals for .^2 and other power operators, as
# there seems to be a bug
normpower(x) = x ./ sqrt.(mean(x.*x,dims=1))
amplify(ratio) = x -> x.*ratio

function ab(tone_len,spacing_len,offset_ratio,repeats,freq,delta,options...)
  @assert 0 <= offset_ratio <= 1
  space = silence(tone_len)

  a_freq=freq
  b_freq=freq*(2^(delta/12))

  a = :without_a in options ? silence(tone_len) : tone(a_freq,tone_len)
  b = :without_b in options ? silence(tone_len) : tone(b_freq,tone_len)
  # a = :ramp in options ? ramp(a,25ms) : a
  # b = :ramp in options ? ramp(b,25ms) : b
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

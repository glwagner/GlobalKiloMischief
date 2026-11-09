using JLD2
for f in ("barotropic_tides_0.36deg_dt10min_convolution_prem_wavedrag3.jld2", "barotropic_tides_0.36deg_dt10min_convolution_paper_wavedrag3.jld2", "barotropic_tides_0.36deg_dt10min_convolution_prem.jld2")
  jldopen(f, "r") do file
    its = sort(parse.(Int, keys(file["timeseries/t"])))
    for it in its[[1:24:end; end]]
      η = file["timeseries/η/$it"][:, :, 1]
      m, I = findmax(abs, η)
      lo = sort(vec(abs.(η)))
      println(f[30:end], " t=", round(file["timeseries/t/$it"]/86400, digits=2), " max ", round(m, digits=2), " at ", Tuple(I), " mean η ", round(sum(η)/length(η), digits=4), " p99 ", round(lo[round(Int, .99length(lo))], digits=2))
    end
  end
end

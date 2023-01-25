import Plots
using LaTeXStrings

Plots.resetfontsizes()
Plots.scalefontsizes(1.5)
Plots.default(titlefont = ("computer modern"),
              legendfont = ("computer modern"),
              linewidth = 2)

ωgrid = 0:.01:2
μgrid = [.1, .6, .8, .9, .99]

ω = [ω for ω in ωgrid, μ in μgrid]
μ = [Complex(μ) for ω in ωgrid, μ in μgrid]

Δ = ω.^4 .* μ.^4 - 4(ω .- 1) .* ω.^2 .* μ.^2
λ₁ = (ω.^2  .* μ.^2 - 2(ω .- 1) + sqrt.(Δ))/2
λ₂ = (ω.^2  .* μ.^2 - 2(ω .- 1) - sqrt.(Δ))/2

Plots.plot(size=(900,450), title=L"|\lambda_+|",
           legend=:bottomright, bottom_margin=3Plots.mm,
           xlabel=L"\omega")
for (i, μ) in enumerate(μgrid)
    # Plots.plot!(ωgrid, abs.(λ₁)[:, i], label="μ=$μ")
    Plots.plot!(ωgrid, abs.(λ₁)[:, i], label="μ=$μ")
end
Plots.plot!()

Plots.savefig("figures/optimal_omega.pdf")

using Plots
using Zygote
using LaTeXStrings

Plots.default(
    # fontfamily="Computer Modern",
    background_color=colorant"#f9f3f1",
    foreground_color=:black,
    linewidth=2,
    xlabelfontsize=18,
    ylabelfontsize=18,
    legendfontsize=18,
    xtickfontsize=18,
    ytickfontsize=18,
    # fontsize=18,
    label=nothing,
)

s = 1
α(x) = .6
F(x) = x - tanh(x - 1)

function animate_fixed_point(F)
    anim = Animation()
    plot(F; label="F(x)", xlims=(-2, 8))
    plot!(y -> y, size=(900, 600), label="x", xlims=(-2, 8))
    x = 7
    for i in 1:10
        newx = F(x)
        println(newx)
        for j in 1:3
            scatter!([x], [x])
            j == 1 && @goto doframe
            plot!([x; x], [x; F(x)])
            scatter!([x; x], [x; F(x)])
            j == 2 && @goto doframe
            plot!([x; F(x)], [F(x); F(x)])
            scatter!([F(x)], [F(x)])
            @label doframe
            frame(anim)
        end
        x = newx
    end
    gif(anim, "fixed_point.gif", fps=2)
end
animate_fixed_point(F)



function animate_chord(f, α, niter=20, fps=2)
    anim = Animation()
    x = 4
    for i in 1:niter
        newx = x - f(x)/α(x)
        println(newx)
        for j in 1:4
            kwargs = (xlims=(-3, 5), ylims=(-1, 1))
            plot(f; kwargs...)
            plot!(y -> 0, size=(900, 600), legend=:none)
            scatter!([x], [f(x)], label=L"x_k")
            j == 1 && @goto doframe
            plot!(y -> f(x) + α(x)*(y - x))
            j == 2 && @goto doframe
            scatter!([newx], [0])
            j == 3 && @goto doframe
            scatter!([newx], [f(newx)])
            j == 4 && @goto doframe
            @label doframe
            frame(anim)
        end
        x = newx
    end
    gif(anim, "newton_raphson.gif", fps=2)
end


s = 1
f(x) = tanh(.3(x - s))
α(x) = 1
animate_chord(f, f', 5)

function animate_bisection(f, a, b, ε = 1e-3)
    anim = Animation()
    i = 0
    while abs(b - a) ≥ ε
        plot(f; xlims=(-.5, 2.4), size=(900, 600))
        plot!([a; b], [0; 0])
        scatter!([a; b], [0; 0])
        annotate!(-.25, 10, text(L"i = " * "$i", :left))
        annotate!(-.25, 9, text(L"|b_i-a_i| = " * "$(abs(b - a))", :left))
        frame(anim)
        x = (a + b) / 2
        scatter!([x], [0])
        a, b = f(a) * f(x) ≤ 0 ? [a, x] : [x, b]
        frame(anim)
        i += 1
    end
    gif(anim, "bisection.gif", fps=2)
end
animate_bisection(x -> x^3 - 2, 0, 2)

function chord(f, α, x; ε = 1e-12, maxiter = 10^6)
    for i in 1:maxiter
        x = x - f(x) / α
        if abs(f(x)) < ε
            return x, i
        end
    end
    error("Failed to converge")
end

using CairoMakie

function f(x, x₀, σ)
    exp(-((x - x₀) / σ)^2)
end

# Boundary condition types
struct Dirichlet end
struct Neumann end

phase(::Dirichlet, cond, val) = cond ? val : -val
phase(::Neumann, cond, val) = val

function evaluate(func, x, L, bc, args...)
    y = mod(x, 2L)
    cond = y < L
    y = cond ? y : 2L - y
    val = func(y, args...)
    phase(bc, cond, val)
end

# Full D'Alembert solution
function sol(x, t, c, L, bc, args...)
    a = x - c * t
    b = x + c * t

    f1 = evaluate(f, a, L, bc, args...)
    f2 = evaluate(f, b, L, bc, args...)
    (f1 + f2) / 2
end


N = 256
L = 1.0
Δx = L / N
xs = (0:N) * Δx
c = 1.0

t = Observable(0.0)

ts = LinRange(0, 1, 256)
figtitle = @lift("t=" * string(round($t, digits=2)))

u_dir = @lift([sol(x, $t, 1, 1, Dirichlet(), 0.5, 0.05) for x in xs])
u_neu = @lift([sol(x, $t, 1, 1, Neumann(), 0.5, 0.05) for x in xs])

fig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=(500, 450), fontsize=24)
    ax1 = Axis(fig[2, 1:3], xlabel=L"x", ylabel=L"u(x, t)", title=L"\text{Dirichlet B.C.: } u(0, t) = u(L, t) = 0")
    ax2 = Axis(fig[3, 1:3], xlabel=L"x", ylabel=L"u(x, t)", title=L"\text{Neumann B.C.: } \partial_x u(0, t) = \partial_x u(L, t) = 0")
    ylims!(ax1, -1.1, 1.1)
    ylims!(ax2, -1.1, 1.1)
    lines!(ax1, xs, u_dir, color=:blue, linewidth=4)
    lines!(ax2, xs, u_neu, color=:red, linewidth=4)
    Label(fig[1, :], figtitle)
    fig
end
##
t[] = 0.55
fig
##

##
record(fig, "test.mp4", ts; framerate=32) do s
    t[] = s
end
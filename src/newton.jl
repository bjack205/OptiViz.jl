
function newton(f, gradf, hessf, x0;
        max_iters = 10,
        tol = 1e-6,
    )
    x = copy(x0)
    for iter = 1:max_iters
        # Compute search direction
        J = f(x)
        ∇f = gradf(x)
        ∇²f = hessf(x)
        δx = -∇²f\∇f


        # Line Search
        ϕ(α) = f(x + α*δx)
        ϕ′(α) = gradf(x + α*δx)'δx
        α = linesearch(ϕ, ϕ′)
        if α != α*0 
            x += α * δx
        end

        # Calculate residual
        res = norm(gradf(x))
        println("Iteration $iter: J = $(f(x)), dJ = $(J - f(x)), α=$α, res=$res")

        if res < tol 
            break
        end

    end
end

function linesearch(ϕ, ϕ′;
        max_iters = 10,
        scaling = 0.5,
        c1 = 1e-4,
        c2 = 0.9
    )
    α = 0.0
    for i = 1:max_iters
        decrease = ϕ(α) ≤ ϕ(0.0) + c1*α*ϕ′(0.0)
        curvature = ϕ′(α) ≥ c2*ϕ′(0.0)
        if decrease && curvature
            return α
        else
            α *= scaling 
        end
    end
    @warn "Line search max iterations"
    α *= 0
    return α
end
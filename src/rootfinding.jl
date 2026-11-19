export quasinewton_rootfinding,bisection_rootfinding

function quasinewton_rootfinding(f,x0, p; tol=1e-4,maxiter=1e4,ϵ=1e-4)
	iter = 0
	x = Inf
    retcode = :Success

	for i in 1:maxiter
		x = x0 - f(x0,p)/df(f,x0,p,ϵ)
		
		abs(x-x0) < tol && return x,retcode
		
		x0 = x
		iter += 1
	end
			
	@warn "Quasi-Newton method did not converge within the specified number of iterations for p = $(p)!"
    retcode = :NotConverged
	return x,retcode
end

function bisection_rootfinding(f,x0,p;Δ = 1.0, a = x0 - Δ, b = x0 + Δ, tol=1e-4, maxiter=1e4)
    retcode = :Success

    if f(a,p) * f(b,p) >= 0
        @warn "Aborting at p = $(p), bisection rootfinding stopped: function values in at a and b must have opposite signs."
        retcode = :NotConverged
        return x0,retcode
    end

    for i in 1:maxiter
        c = (a + b) / 2  # Calculate the midpoint of the interval
        fc = f(c,p)        # Evaluate the function at the midpoint

        if abs(fc) < tol || (b - a) / 2 < tol
            return c,retcode
        end

        if f(a,p) * fc < 0
            b = c
        else
            a = c
        end
    end

    @warn "Bisection method did not converge within the specified number of iterations p = $(p)!"
    retcode = :NotConverged
	return x,retcode
end


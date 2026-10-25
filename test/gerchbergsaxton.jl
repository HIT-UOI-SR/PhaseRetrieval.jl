@testset "Gerchberg–Saxton" begin
    rng = MersenneTwister(23)
    truth = randn(rng, ComplexF64, 7, 6)
    ain, aout = abs.(truth), abs.(fft(truth))
    init = randn(rng, ComplexF64, size(truth))
    saved = copy(init)
    expected = amplitude(ain, ifft(amplitude(aout, fft(amplitude(ain, init)))))
    x, y = gerchbergsaxton(ain, aout; init=init, maxiter=1, tol=0)
    @test x ≈ expected
    @test y ≈ fft(x)
    @test abs.(x) ≈ ain
    @test init == saved

    x0, y0 = gerchbergsaxton(ain, aout; init=init, maxiter=0)
    @test x0 ≈ amplitude(ain, init)
    @test y0 ≈ fft(x0)
    x, y = gerchbergsaxton(ain, aout; init=init, maxiter=300, tol=0)
    @test msd(abs.(y), aout) < 0.02msd(abs.(y0), aout)
    # Exact solutions remain fixed, including masked zeros and global phase.
    truth[1:2, :] .= 0
    x, y = gerchbergsaxton(abs.(truth), abs.(fft(truth)); init=truth .* cis(0.8), maxiter=5, tol=0)
    @test x ≈ truth .* cis(0.8)
    @test y ≈ fft(x)

    # A custom unitary transform need not be Fourier, square, or shape preserving.
    U = Matrix(qr(randn(rng, ComplexF64, 8, 8)).Q)
    target = randn(rng, ComplexF64, 2, 4)
    f = x -> U * vec(x)
    b = y -> reshape(U' * y, 2, 4)
    x, y = gerchbergsaxton(abs.(target), abs.(f(target)); projector=(f,b), init=target, tol=0, maxiter=2)
    @test x ≈ target
    @test y ≈ f(x)

    # Stopping checks the predicted amplitude rather than the imposed one.
    calls = Ref(0)
    counted = x -> (calls[] += 1; fft(x))
    gerchbergsaxton(ain, aout; init=init, projector=(counted,ifft), maxiter=8, tol=0)
    @test calls[] == 9
    calls[] = 0
    gerchbergsaxton(ain, aout; init=init, projector=(counted,ifft), maxiter=8, tol=1e100)
    @test calls[] == 1
end

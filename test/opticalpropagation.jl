@testset "OpticalPropagation 0.4" begin
    # An offset, rectangular, odd/even grid exercises the physical-grid metadata.
    data = [exp(-3 * (x^2 + y^2)) * cis(0.8x^2 - 0.6y)
            for x in range(-1, 1; length=9), y in range(-1, 1; length=8)]
    truth = MonoLightField2D(data; wavelength=532u"nm",
                            size=(0.8u"mm", 0.6u"mm"),
                            center=(0.2u"mm", -0.1u"mm"))
    distance = 20u"mm"

    for method in (fresnel2, angularspectrum)
        @testset "$method" begin
            forward, backward = method(distance), method(-distance)
            measured = forward(truth)

            # A solution must remain a fixed point of the GS iterations.
            r1, r2 = gerchbergsaxton(abs.(truth), abs.(measured);
                                    projector=(forward, backward), init=copy(truth),
                                    maxiter=3, tol=0)
            @test r1 ≈ truth
            @test r2 ≈ measured
            for result in (r1, r2)
                @test result isa MonoLightField2D
                @test physicalsize(result) == physicalsize(truth)
                @test gridcenter(result) == gridcenter(truth)
                @test samplepitch(result) == samplepitch(truth)
                @test result.wavelength == truth.wavelength
            end

            # Retrieve from an incorrect initial phase, using ordinary matrices.
            matrixforward = x -> parent(forward(MonoLightField2D(truth; data=x)))
            matrixbackward = x -> parent(backward(MonoLightField2D(truth; data=x)))
            ain, aout = abs.(parent(truth)), abs.(parent(measured))
            init = complex.(ain)
            initial_error = msd(abs.(matrixforward(init)), aout)
            r1, r2 = gerchbergsaxton(ain, aout;
                                    projector=(matrixforward, matrixbackward),
                                    init=copy(init), maxiter=30, tol=0)
            @test abs.(r1) ≈ ain
            @test r2 ≈ matrixforward(r1)
            @test msd(abs.(matrixforward(r1)), aout) < initial_error
        end
    end
end

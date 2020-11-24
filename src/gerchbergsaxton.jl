"""
    gerchbergsaxton(Ain, Aout, kw...)

reconstruct the complex signals from their amplitude values based on the Gerchberg-Saxton algorithm.

## References

1. R. W. Gerchberg and W. O. Saxton, "A practical algorithm for the determination of the phase from image and diffraction plane pictures," Optik 35, 237 (1972)
"""
function gerchbergsaxton(Ain, Aout; 
        maxiter=100,
        projector=(fft, ifft),
        init=@. exp(2π * im * $(rand(Float64, size(Ain))))
    )
    Cin = init
    Cout = similar(Aout)
    for i in 1:maxiter
        @. Cin = Ain * exp(2π * im * angle(Cin))
        Cout = projector[1](Cin)
        @. Cout = Aout * exp(2π * im * angle(Cout))
        Cin = projector[2](Cout)
    end
    (Cin, Cout)
end

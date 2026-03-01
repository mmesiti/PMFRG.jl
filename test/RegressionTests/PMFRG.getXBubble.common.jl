function compare_arguments_post(Xexp, X)
    result = true
    for field in fieldnames(typeof(X))
        val = getfield(X, field)
        expval = getfield(Xexp, field)
        # This might be too strict
        @assert typeof(expval) == Array{Float64,4}
        if val != expval
            abs_diff = abs.(val .- expval)
            total_diff = sum(abs_diff)
            if total_diff > 1e-14
                different_val_places = val .!= expval

                print("Test failed: $field differs significantly: ")
                print("Absolute Difference: ", total_diff, ", ")
                print("Maximum Difference: ", maximum(abs_diff), ", ")

                println("First computed values: ")
                println(first(val[different_val_places], 5))
                println("First computed values - reinterpreted: ")
                println(reinterpret(UInt64, first(val[different_val_places], 5)))

                println("First expected values: ")
                println(first(expval[different_val_places], 5))
                println("First expected values - reinterpreted: ")
                println(reinterpret(UInt64, first(expval[different_val_places], 5)))
                result = false
            end
        end
    end
    result
end

using SpinFRGLattices.SquareLattice
function generate_test_params()
    System = getSquareLattice(4, [1.0, 0.1])
    Params(System, OneLoop(), T = 0.5, N = 5, accuracy = 1e-3)
end




using HDF5
using H5Zblosc
function h5deserialize(f::HDF5.File, root::String, id::Int)
    (
        X = PMFRG.BubbleType(
            read(f["$root/$id/X/a"]),
            read(f["$root/$id/X/b"]),
            read(f["$root/$id/X/c"]),
            read(f["$root/$id/X/Ta"]),
            read(f["$root/$id/X/Tb"]),
            read(f["$root/$id/X/Tc"]),
            read(f["$root/$id/X/Td"]),
        ),
        State = PMFRG.StateType(
            read(f["$root/$id/State/f_int"]),
            read(f["$root/$id/State/γ"]), # This might need change
            PMFRG.HeisenbergVertex(
                read(f["$root/$id/State/Γ/a"]),# This might need change
                read(f["$root/$id/State/Γ/b"]),# This might need change
                read(f["$root/$id/State/Γ/c"]),# This might need change
            ),
        ),
        Deriv = PMFRG.StateType(
            read(f["$root/$id/Deriv/f_int"]),
            read(f["$root/$id/Deriv/γ"]), # This might need change
            PMFRG.HeisenbergVertex(
                read(f["$root/$id/Deriv/Γ/a"]),# This might need change
                read(f["$root/$id/Deriv/Γ/b"]),# This might need change
                read(f["$root/$id/Deriv/Γ/c"]),# This might need change
            ),
        ),
        Lam = read(f["$root/$id/Lam"]),
    )
end

"""Function that can be used to store X,State,Deriv and Lam
(the non-trivially-constructible arguments of getXBubble!)
in a HDF5 file, so that h5deserialize can read them properly"""
function h5serialize!(
    f::HDF5.File,
    root::String,
    X::PMFRG.BubbleType,
    S::PMFRG.StateType,
    D::PMFRG.StateType,
    Lam::Float64,
    id::Int,
)
    # X
    f["$root/$id/X/a", blosc = 9] = X.a
    f["$root/$id/X/b", blosc = 9] = X.b
    f["$root/$id/X/c", blosc = 9] = X.c
    f["$root/$id/X/Ta", blosc = 9] = X.Ta
    f["$root/$id/X/Tb", blosc = 9] = X.Tb
    f["$root/$id/X/Tc", blosc = 9] = X.Tc
    f["$root/$id/X/Td", blosc = 9] = X.Td

    # S
    f["$root/$id/State/f_int", blosc = 9] = S.f_int
    f["$root/$id/State/iSigma", blosc = 9] = S.iSigma
    f["$root/$id/State/Γ/a", blosc = 9] = S.Gamma.a
    f["$root/$id/State/Γ/b", blosc = 9] = S.Gamma.b
    f["$root/$id/State/Γ/c", blosc = 9] = S.Gamma.c

    # D
    f["$root/$id/Deriv/f_int", blosc = 9] = D.f_int
    f["$root/$id/Deriv/iSigma", blosc = 9] = D.iSigma # This might need change
    f["$root/$id/Deriv/Γ/a", blosc = 9] = D.Gamma.a
    f["$root/$id/Deriv/Γ/b", blosc = 9] = D.Gamma.b
    f["$root/$id/Deriv/Γ/c", blosc = 9] = D.Gamma.c


    # Lam
    f["$root/$id/Lam"] = Lam
end


"""This function has only been used
to convert the existing serialized data
into HDF5.
It is kept here only for reference"""
function convert_jld2_to_hdf5(data, filename)

    println("Saving to hdf5")
    f = h5open(filename, "w")
    try
        ncases = length(data["return_value"])
        f["Ncases"] = ncases
        for i = 1:ncases
            workspace, lam, _ = (data["arguments"])[i]
            workspace_post_exp, _, _ = (data["arguments_post"])[i]

            let X, State, Deriv
                (; X, State, Deriv) = workspace
                h5serialize!(f, "arguments", X, State, Deriv, lam, i)
            end


            let X, State, Deriv
                (; X, State, Deriv) = workspace_post_exp
                h5serialize!(f, "arguments_post", X, State, Deriv, lam, i)
            end


        end
    finally
        close(f)
    end

end

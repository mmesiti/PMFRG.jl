"""Exectutes nontrivial symmetry in flow equations for Heisenberg dimer. Local and nonlocal b vertex are equal: Gammab_11 = Gammab_12!"""
function test_runFRG(Method = OneLoop(), System = getPolymer(2))
    Par = BenchmarkingParams(Method, System)
    test_runFRG(Par)
end

function test_runFRG(Par::PMFRG.PMFRGParams; kwargs...)
    tempFolder = "temp_PMFRG_test"
    mainFile = joinpath(tempFolder, "temp_main.h5")
    CheckPoints = joinpath(tempFolder, "Checkpoints.h5")

    SolP, ObsPt =
        SolveFRG(Par, MainFile = mainFile, CheckpointDirectory = CheckPoints; kwargs...)

    println("cleaning up... deleting ", mainFile, " and ", CheckPoints)
    rm(tempFolder, recursive = true)
    iSigma = SolP.u[end].x[2]
    Gammaa = SolP.u[end].x[3]
    Gammab = SolP.u[end].x[4]
    Gammac = SolP.u[end].x[5]
    return iSigma, Gammaa, Gammab, Gammac, ObsPt.saveval[end], ObsPt.t, Par
end


function test_DimerFRG(Method = OneLoop(); Obsacc = 1e-14, kwargs...)
    iSigma, Gammaa, Gammab, Gammac, Obs, t, Par = test_runFRG(Method, getPolymer(2))
    test_FRGResults(iSigma, Gammaa, Gammab, Gammac, Obs, t, Par; Obsacc = Obsacc, kwargs...)

    println(
        "Testing whether local and nonlocal b vertex are equal on dimer Gammab_11 = Gammab_12",
    )
    test_Gammab_Dimer(Gammab; kwargs...)

end

function test_SquagomeFRG(Method = OneLoop(); Obsacc = 1e-14, kwargs...)
    SysFunc = SquareKagome.getMirrorSquareKagome
    iSigma, Gammaa, Gammab, Gammac, Obs, t, Par = test_runFRG(Method, SysFunc(4, 1, 0.2))
    test_FRGResults(
        iSigma,
        Gammaa,
        Gammab,
        Gammac,
        Obs,
        t,
        Par,
        (Method, SysFunc);
        Obsacc = Obsacc,
        kwargs...,
    )

end

function test_FRGResults(
    iSigma,
    Gammaa,
    Gammab,
    Gammac,
    Obs,
    t,
    Par,
    Method = PMFRG.getPMFRGMethod(Par);
    Obsacc = 1e-14,
    kwargs...,
)
    @testset "Λ saved values" begin
        @test t[end] ≈ Par.NumericalParams.Lam_min
    end

    test_nonzero(iSigma, Gammaa, Gammab, Gammac, Par.System.couplings)

    @testset verbose = true "Testing frequency symmetries" begin
        println("Testing onsite Gammaa_ii vertex")
        test_Gammaa_onsite(Gammaa; kwargs...)

        println("Testing t ↔ u symmetries")
        test_tu_symmetries(Gammaa, Gammab, Gammac; kwargs...)
    end

    test_Observables(Method, Obs, Obsacc = Obsacc)

    return
end

function test_runDimerParquet()
    ParquetLambda = 0.0
    tempFolder = "temp_PMFRG_test"

    MainFile = joinpath(tempFolder, "temp_main.h5")

    CheckpointDirectory = joinpath(tempFolder, "Checkpoints.h5")


    Par = BenchmarkingParams(Parquet())
    Sol, Obs = SolveParquet(Par, ParquetLambda; MainFile, CheckpointDirectory)
    println("cleaning up... deleting ", MainFile, " and ", CheckpointDirectory)
    rm(tempFolder, recursive = true)
    return Sol, Obs, Par
end

function test_DimerParquet(; kwargs...)
    Sol, Obs, Par = test_runDimerParquet()
    Gammaa, Gammab, Gammac = Sol.State.Gamma.a, Sol.State.Gamma.b, Sol.State.Gamma.c
    test_BareBubbles(Sol; kwargs...)
    test_BubbleSymmetries(Sol; kwargs...)
    test_FRGResults(
        Sol.State.iSigma,
        Gammaa,
        Gammab,
        Gammac,
        Obs,
        Par.NumericalParams.Lam_min,
        Par;
        kwargs...,
    )
end

function test_Observables(Method, Obs; Obsacc = 1e-14)
    obs_ex = example_Obs(Method)
    function test(ObsName)
        O1 = getproperty(Obs, ObsName)
        O2 = getproperty(obs_ex, ObsName)
        for i in eachindex(O2)
            @test O1[i] ≈ O2[i] atol = Obsacc
        end
    end

    println("Observables: ", Obs)
    @testset "Testing Susceptibility" begin
        test(:Chi)
    end

    @testset "Testing iSigma" begin
        test(:iSigma)
    end

    @testset "Testing max(Gamma)" begin
        @testset "Gammaa" begin
            test(:MaxVa)
        end
        @testset "Gammab" begin
            test(:MaxVb)
        end
        @testset "Gammac" begin
            test(:MaxVc)
        end
    end
end

test_Observables(Method::Parquet, Obs; Obsacc = 1e-14) =
    println("Observables check not implemented for parquet")

function test_nonzero(iSigma, Gammaa, Gammab, Gammac, couplings)
    @testset "non-trivial Vertices" begin
        @testset "iSigma" begin
            @test maximum(abs, iSigma) > 0.0
        end
        @testset "Gammaa" begin
            @test maximum(abs, Gammaa) > 0.0
        end
        @testset "Gammab" begin
            @test maximum(abs, Gammab) > 0.0
        end
        @testset "Gammac" begin
            @test maximum(abs, Gammac) > maximum(abs, couplings)
        end
    end
end


function test_tu_symmetries(Gammaa, Gammab, Gammac; kwargs...)
    test_tu_symmetry_ab(Gammaa, "Gammaa"; kwargs...)
    test_tu_symmetry_ab(Gammab, "Gammab"; kwargs...)
    test_tu_symmetry_c(Gammaa, Gammab, Gammac; kwargs...)
end

function test_Gammab_Dimer(Gammab::AbstractArray; tol = 1e-14)
    vb1 = @view Gammab[1, :, :, :]
    vb2 = @view Gammab[2, :, :, :]
    @testset "Dimer: Gammab_11 == Gammab_12" begin
        Failures = 0
        for i in eachindex(vb1, vb2)
            if !isapprox(vb1[i], vb2[i], atol = tol)
                Failures += 1
            end
        end
        @test Failures == 0
    end
end

function test_Gammaa_onsite(Gammaa::AbstractArray, OnsiteBonds = [1]; tol = 1e-14)
    @testset "Gammaa_ii(stu) == -Gammaa_ii(uts)" begin
        Failures = 0
        for Rij in OnsiteBonds,
            s in axes(Gammaa, 2),
            t in axes(Gammaa, 3),
            u in axes(Gammaa, 4)

            if !isapprox(Gammaa[Rij, s, t, u], -Gammaa[Rij, u, t, s], atol = tol)
                Failures += 1
            end
        end
        @test Failures == 0
    end

    @testset "Gammaa_ii(stu) == Gammaa_ii(tus)" begin
        Failures = 0
        for Rij in OnsiteBonds,
            s in axes(Gammaa, 2),
            t in axes(Gammaa, 3),
            u in axes(Gammaa, 4)

            if !isapprox(Gammaa[Rij, s, t, u], +Gammaa[Rij, t, u, s], atol = tol)
                Failures += 1
            end
        end
        @test Failures == 0
    end
end

function test_tu_symmetry_ab(Gamma::AbstractArray, Name; tol = 1e-14)

    @testset "$(Name)_ij(stu) == -$(Name)_ij(sut)" begin
        Failures = 0
        for Rij in axes(Gamma, 1),
            s in axes(Gamma, 2),
            t in axes(Gamma, 3),
            u in axes(Gamma, 4)

            if !isapprox(Gamma[Rij, s, t, u], -Gamma[Rij, s, u, t], atol = tol)
                Failures += 1
            end
        end
        @test Failures == 0

    end
end

function test_tu_symmetry_c(
    Gammaa::AbstractArray,
    Gammab::AbstractArray,
    Gammac::AbstractArray,
    V = "Gamma";
    tol = 1e-14,
)
    Failures = 0
    @testset "$(V)c_ij(stu) == (-$(V)a_ij + $(V)b_ij + $(V)c_ij)(sut)" begin
        for Rij in axes(Gammac, 1),
            s in axes(Gammac, 2),
            t in axes(Gammac, 3),
            u in axes(Gammac, 4)

            if !isapprox(
                Gammac[Rij, s, t, u],
                (-Gammaa[Rij, s, u, t] + Gammab[Rij, s, u, t] + Gammac[Rij, s, u, t]),
                atol = tol,
            )
                Failures += 1
            end
        end
        @test Failures == 0
    end

end

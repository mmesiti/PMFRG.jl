import DiffEqCallbacks: SavedValues

function test_saving(
    FileDirectory,
    Method = OneLoop(),
    GeometryGenerator = SquareKagome.getSquareKagome,
)
    FileDirectory = PMFRGCore.UniqueDirName(FileDirectory)
    mkpath(FileDirectory)
    Par = BenchmarkingParams(Method, GeometryGenerator(4))
    State = PMFRGCore.InitializeState(Par)
    Lam = Par.NumericalParams.Lam_min + 0.01
    saved_values = SavedValues(Float64, PMFRGCore.Observables)

    obs1 = PMFRGCore.getObservables(
        PMFRGCore.Observables,
        State,
        Par.NumericalParams.Lam_max,
        Par,
    )
    obs2 = PMFRGCore.getObservables(PMFRGCore.Observables, State, Lam, Par)

    push!(saved_values.t, Par.NumericalParams.Lam_max)
    push!(saved_values.saveval, obs1)

    push!(saved_values.t, Lam)
    push!(saved_values.saveval, obs2)

    Filename = PMFRGCore.saveCurrentState(FileDirectory, State, saved_values, Lam, Par)
    return Filename
end

function test_IO(Method = OneLoop(), GeometryGenerator = SquareKagome.getSquareKagome)
    @testset "FileIO" begin
        tempFolder = joinpath("temp_PMFRG_test")
        CheckpointDir = joinpath("temp_PMFRG_test", "Checkpoints")
        Filename = test_saving(CheckpointDir, Method, GeometryGenerator)

        Par = PMFRGCore.readParams(Filename, GeometryGenerator)
        @testset "Testing loopOrder from Params" begin
            @test PMFRGCore.getPMFRGMethod(Par) == Method
        end
        rm(tempFolder, recursive = true)
        println("cleaned up... deleted ", tempFolder)
    end
end

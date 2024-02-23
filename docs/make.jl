using Documenter, PMFRG
#push!(LOAD_PATH,"../src/")

makedocs(sitename="PMFRG.jl", repo=Remotes.GitHub("mmesiti","PMFRG.jl"))
#deploydocs(;repo="git@gitlabph.physik.fu-berlin.de:niggeni/PMFRG.jl",)
deploydocs(;repo="github.com/mmesiti/PMFRG.jl.git",)
# deploydocs(;repo="github.com/NilsNiggemann/PMFRG.jl",)

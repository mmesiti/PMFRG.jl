using Documenter, PMFRG
#push!(LOAD_PATH,"../src/")

makedocs(sitename="PMFRG.jl",
         repo=Remotes.GitHub("mmesiti","PMFRG.jl"),
         format=Documenter.HTML(edit_link=:commit))
deploydocs(repo="github.com/mmesiti/PMFRG.jl",devbranch="docs")

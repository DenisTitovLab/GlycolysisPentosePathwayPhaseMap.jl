using Test, GlycolysisPentosePathwayPhaseMap
const G = GlycolysisPentosePathwayPhaseMap

@testset "update_provenance! creates, replaces one table, preserves the other" begin
    path = joinpath(mktempdir(), "provenance.toml")
    G.update_provenance!(path, "grids", Dict("ppp_commit" => "abc", "atp_grid_rows" => 3))
    @test G.read_provenance(path)["grids"]["ppp_commit"] == "abc"

    G.update_provenance!(path, "literature", Dict("ppp_experiments_commit" => "def"))
    p = G.read_provenance(path)
    @test p["grids"] == Dict("ppp_commit" => "abc", "atp_grid_rows" => 3)
    @test p["literature"]["ppp_experiments_commit"] == "def"

    G.update_provenance!(path, "grids", Dict("ppp_commit" => "xyz"))
    p = G.read_provenance(path)
    @test p["grids"] == Dict("ppp_commit" => "xyz")          # old keys dropped
    @test p["literature"] == Dict("ppp_experiments_commit" => "def")   # untouched
end

@testset "git_state" begin
    @test G.git_state(mktempdir()) === nothing               # not a git work tree
    repo = mktempdir()
    run(pipeline(`git -C $repo init -q`; stderr = devnull))
    write(joinpath(repo, "a.txt"), "1")
    run(`git -C $repo -c user.name=t -c user.email=t@t add a.txt`)
    run(`git -C $repo -c user.name=t -c user.email=t@t commit -q -m init`)
    s = G.git_state(repo)
    @test length(s.commit) == 40
    @test !s.dirty
    write(joinpath(repo, "b.txt"), "2")
    @test G.git_state(repo).dirty                            # untracked file -> dirty
    @test !G.git_state(repo; paths = ["a.txt"]).dirty        # scoped to clean path
    @test G.git_state(joinpath(repo, "sub_does_not_exist")) === nothing
end

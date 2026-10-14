# provenance.jl — where the shipped data came from (data/provenance.toml)
#
# Two writers share one file, each owning one table: `regenerate_grid` (the recompute
# extension) writes [grids]; the PPP_Experiments literature export writes [literature].
# update_provenance! read-modify-writes only the named table so neither clobbers the other.

"Path of the shipped provenance file."
provenance_path() = joinpath(DATA_DIR, "provenance.toml")

"Parse a provenance TOML file into a `Dict{String,Any}`."
read_provenance(path = provenance_path()) = TOML.parsefile(path)

"""
    update_provenance!(path, section, table) -> Dict{String,Any}

Replace table `section` of the provenance TOML at `path` with `table` (creating the file if
absent), leaving every other table untouched. Keys are written sorted. Returns the full
parsed content as written.
"""
function update_provenance!(path::AbstractString, section::AbstractString, table::AbstractDict)
    prov = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
    prov[String(section)] = Dict{String,Any}(String(k) => v for (k, v) in table)
    mkpath(dirname(path))
    open(io -> TOML.print(io, prov; sorted = true), path, "w")
    return prov
end

"""
    git_state(dir; paths = String[]) -> Union{Nothing, NamedTuple}

`(; commit, branch, dirty)` of the git work tree whose top level is `dir`, or `nothing` when
`dir` is not a work-tree top level (e.g. a registry-installed package). `dirty` counts
untracked files too and is restricted to `paths` (relative to `dir`) when given.
"""
function git_state(dir::AbstractString; paths = String[])
    isdir(dir) || return nothing
    git(args...) = readchomp(pipeline(`git -C $dir $args`; stderr = devnull))
    top = try git("rev-parse", "--show-toplevel") catch; return nothing end
    realpath(top) == realpath(dir) || return nothing
    commit = git("rev-parse", "HEAD")
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    dirty  = !isempty(git("status", "--porcelain", "--", paths...))
    return (; commit, branch, dirty)
end

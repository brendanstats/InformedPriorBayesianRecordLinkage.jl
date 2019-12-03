function findnz_dim(A::Array{T, 3}, dim::Integer) where T <: Real
    sparseA = sparse(A[:, :, dim])
    return hcat(findnz(sparseA)..., fill(T(dim), nnz(sparseA)))
end

findnz_stagect(A::Array{T, 3}) where T <: Real = mapreduce(ii -> findnz_dim(A, ii), vcat, 1:size(A, 3))

function h5write_BlockStages(filename::String,
                             bstages::BlockStages;
                             groupname::String = "/",
                             mode::String = "w")
    h5open(filename, mode) do writef
        writef[groupname * "/rs2b"] = bstages.rs2b
        writef[groupname * "/cs2b"] = bstages.cs2b
        writef[groupname * "/b2rows"] = bstages.b2rows
        writef[groupname * "/n2cols"] = bstages.b2cols
        writef[groupname * "/b2conditionedb"] = bstages.b2conditionedb
        writef[groupname * "/b2s"] = bstages.b2s
        writef[groupname * "/bnrows"] = bstages.bnrows
        writef[groupname * "/bncols"] = bstages.bncols
        writef[groupname * "/bmindim"] = bstages.bmindim
        writef[groupname * "/bmaxdim"] = bstages.bmaxdim
        writef[groupname * "/smindim"] = bstages.smindim
        writef[groupname * "/nrow"] = bstages.nrow
        writef[groupname * "/ncol"] = bstages.ncol
        writef[groupname * "/nblock"] = bstages.nblock
        writef[groupname * "/nstage"] = bstages.nstage
    end
    nothing
end

function h5read_BlockStages(filename::String, groupname::String = "/")
    return h5open(filename, "r") do readf
        BlockStages(
            readf[groupname * "/rs2b"],
            readf[groupname * "/cs2b"],
            readf[groupname * "/b2rows"],
            readf[groupname * "/n2cols"],
            readf[groupname * "/b2conditionedb"],
            readf[groupname * "/b2s"],
            readf[groupname * "/bnrows"],
            readf[groupname * "/bncols"],
            readf[groupname * "/bmindim"],
            readf[groupname * "/bmaxdim"],
            readf[groupname * "/smindim"],
            readf[groupname * "/nrow"],
            readf[groupname * "/ncol"],
            readf[groupname * "/nblock"],
            readf[groupname * "/nstage"]
        )
    end
end

function h5write_MCMC(filename::String,
                      pairlinkcounts::SparseMatrixCSC,
                      nlinks::Array{<:Integer, 1},
                      pM::Array{<:Real, 2},
                      pU::Array{<:Real, 2},
                      linkmoves::Integer,
                      nsteps::Integer = length(nlinks);
                      groupname::String = "/",
                      mode::String = "w")
    h5open(filename, mode) do writef
        writef[groupname * "/pairlinkcounts"] = hcat(findnz(pairlinkcounts)...)
        writef[groupname * "/nlinks"] = nlinks
        writef[groupname * "/pM"] = pM
        writef[groupname * "/pU"] = pU
        writef[groupname * "/linkmoves"] = linkmoves
        writef[groupname * "/nsteps"] = nsteps
    end
    nothing
end

function h5write_StageMCMC(filename::String,
                           pairlinkcounts::Array{<:Integer, 3},
                           nlinks::Array{<:Integer, 2},
                           pM::Array{<:Real, 2},
                           pU::Array{<:Real, 2},
                           linkmoves::Integer,
                           stagemovese::Integer,
                           nsteps::Integer = size(nlinks, 1);
                           groupname::String = "/",
                           mode::String = "w")
    h5open(filename, mode) do writef
        writef[groupname * "/pairlinkcounts"] = findnz_stagect(pairlinkcounts)
        writef[groupname * "/nlinks"] = nlinks
        writef[groupname * "/pM"] = pM
        writef[groupname * "/pU"] = pU
        writef[groupname * "/linkmoves"] = linkmoves
        writef[groupname * "/stagemoves"] = stagemoves
        writef[groupname * "/nsteps"] = nsteps
    end
    nothing
end

function h5write_StageMCMC(filename::String,
                           pairlinkcounts::Array{<:Integer, 3},
                           nlinks::Array{<:Integer, 2},
                           linkmoves::Integer,
                           stagemoves::Integer,
                           nsteps::Integer = size(nlinks, 1);
                           groupname::String = "/",
                           mode::String = "w")
    h5open(filename, mode) do writef
        writef[groupname * "/pairlinkcounts"] = findnz_stagect(pairlinkcounts)
        writef[groupname * "/nlinks"] = nlinks
        writef[groupname * "/linkmoves"] = linkmoves
        writef[groupname * "/stagemoves"] = stagemoves
        writef[groupname * "/nsteps"] = nsteps
    end
    nothing
end

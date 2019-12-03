struct BlockStages{G <: Integer}
    rs2b::Array{G, 2}
    cs2b::Array{G, 2}
    b2rows::Dict{G, Array{G, 1}}
    b2cols::Dict{G, Array{G, 1}}
    b2conditionedb::Dict{G, Array{G, 1}}
    b2s::Array{G, 1}
    bnrows::Array{G, 1}
    bncols::Array{G, 1}
    bmindim::Array{G, 1}
    bmaxdim::Array{G, 1}
    smindim::Array{G, 1}
    nrow::G
    ncol::G
    nblock::G
    nstage::G
end

"""
rs2b = [[1, 1, 2, 2, 3] [4, 4, 4, 4, 4]]
cs2b = [[1, 2, 3, 3] [4, 4, 4, 4]]
nblock = 4
G = Int64

bstages = BlockStages(rs2b, cs2b)

bstages.rs2b == rs2b
bstages.cs2b == cs2b
bstages.b2rows == Dict(1 => [1, 2], 2 => [3, 4], 3 => [5], 4 => [1, 2, 3, 4, 5])
bstages.b2cols == Dict(1 => [1], 2 => [2], 3 => [3, 4], 4 => [1, 2, 3, 4])
bstages.b2s == [1, 1, 1, 2]
bstages.b2nextb == [4, 4, 4, 0]
bstages.bnrows == [2, 2, 1, 5]
bstages.bncols == [1, 1, 2, 4]
bstages.bmindim == [1, 1, 1, 4]
bstages.nrow == 5
bstages.ncol == 4
bstages.nblock == 4
bstages.nstage == 2

bstages = BlockStages(Int32.(rs2b), Int32.(cs2b))
"""
function BlockStages(rs2b::Array{G, 2}, cs2b::Array{G, 2}, nblock::Integer = maximum(rs2b)) where {G <: Integer}
    nstage = G(size(rs2b, 2))
    nrow = G(size(rs2b, 1))
    ncol = G(size(cs2b, 1))

    b2rows = Dict{G, Array{G, 1}}()
    b2cols = Dict{G, Array{G, 1}}()
    for kk in 1:nblock
        b2rows[kk] = G[]
        b2cols[kk] = G[]
    end

    for ii in one(G):nrow, ss in one(G):nstage
        if !iszero(rs2b[ii, ss])
            push!(b2rows[rs2b[ii, ss]], ii)
        end
    end

    for jj in one(G):ncol, ss in one(G):nstage
        if !iszero(cs2b[jj, ss])
            push!(b2cols[cs2b[jj, ss]], jj)
        end
    end
    
    b2s = ones(G, nblock)
    if nstage > one(G)
        for ss in 2:nstage
            for bb in unique(rs2b[:, ss])
                if !iszero(bb) #allow rows to not be present in some stages
                    b2s[bb] = ss
                end
            end
        end
    end

    b2conditionedb = Dict{G, Array{G, 1}}()
    for bb in one(G):nblock
        if b2s[bb] < nstage
            #assume all blocks have at least one row
            b2conditionedb[bb] = rs2b[b2rows[bb][1], (b2s[bb] + one(G)):nstage]
        else
            b2conditionedb[bb] = G[]
        end
    end

    bnrows = map(bb -> G(length(b2rows[bb])), one(G):nblock)
    bncols = map(bb -> G(length(b2cols[bb])), one(G):nblock)
    bmindim = zeros(G, nblock)
    bmaxdim = zeros(G, nblock)
    smindim = zeros(G, nstage)
    for block in one(G):nblock
        if bnrows[block] >= bncols[block]
            bmaxdim[block] = bnrows[block]
            bmindim[block] = bncols[block]
            smindim[b2s[block]] += bncols[block]
        else
            bmaxdim[block] = bncols[block]
            bmindim[block] = bnrows[block]
            smindim[b2s[block]] += bnrows[block]
        end
    end
    
    return BlockStages(rs2b, cs2b, b2rows, b2cols, b2conditionedb, b2s, bnrows, bncols, bmindim, bmaxdim, smindim, nrow, ncol, nblock, nstage)
end

function BlockStages(nrow::G, ncol::G) where G <: Integer
    return BlockStages(fill(one(G), (nrow, one(G))), fill(one(G), (ncol, one(G))), one(G))
end

function minimumstage_link(row::Integer, col::Integer, bstages::BlockStages)
    return findfirst(map(stage -> iszero(bstages.rs2b[row, stage]) ? false : bstages.rs2b[row, stage] == bstages.cs2b[col, stage], 1:bstages.nstage))
end

#bstages.rs2b[row, :] to get sequence of blocks

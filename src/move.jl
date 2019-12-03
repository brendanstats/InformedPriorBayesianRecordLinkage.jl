"""
Count number of links excluded from each block due to conditioning
"""
function count_blockprevlinks(C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer}
    blockprevlinks = zeros(G, bstages.nblock)
    for block in one(G):bstages.nblock
        if !iszero(bstages.b2nextb[block])
            blockprevlinks[bstages.b2nextb[block]] += C.nlinkblock[block] + blockprevlinks[block]
        end
    end
    return blockprevlinks
end

"""
Count number of links excluded from each stage due to conditioning
"""
function count_stageprevlinks(C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer}
    stageprevlinks = zeros(G, bstages.nstage)
    if bstages.nstage > one(G)
        for stage in one(G):(bstages.nstage - one(G))
            stageprevlinks[stage + one(G)] = stageprevlinks[stage] + C.nlinkstage[stage]
        end
    end
    return stageprevlinks
end

"""
Count number of links excluded from each block and each stage due to conditioning
"""
function count_prevlinks(C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer}
    return count_blockprevlinks(C, bstages), count_stageprevlinks(C, bstages)
end

function loglik_add(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf))  where {T <: AbstractFloat}
    loglikratio = get_loglik(row, col, compsum, loglikRatios, loglikMissing)
    return loglikratio
end

function loglik_remove(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf))  where {T <: AbstractFloat}
    loglikratio = -get_loglik(row, col, compsum, loglikRatios, loglikMissing)
    return loglikratio
end

function loglik_rowswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf))  where {T <: AbstractFloat}
    loglikratio = get_loglik(row, col, compsum, loglikRatios, loglikMissing) - get_loglik(C.col2row[col], col, compsum, loglikRatios, loglikMissing)
    return loglikratio
end

function loglik_colswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf))  where {T <: AbstractFloat}
    loglikratio = get_loglik(row, col, compsum, loglikRatios) - get_loglik(row, C.row2col[row], compsum, loglikRatios, loglikMissing)
    return loglikratio
end

function loglik_doubleswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf))  where {T <: AbstractFloat}
    rowalt = C.col2row[col]
    colalt = C.row2col[row]
    loglikratio = get_loglik(row, col, compsum, loglikRatios, loglikMissing)
    loglikratio += get_loglik(rowalt, colalt, compsum, loglikRatios, loglikMissing)
    loglikratio -= get_loglik(row, colalt, compsum, loglikRatios, loglikMissing)
    loglikratio -= get_loglik(rowalt, col, compsum, loglikRatios, loglikMissing)
    return loglikratio
end

function counts_add(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2})
    countsdelta = get_counts(row, col, compsum, obsidxCounts)
    return countsdelta
end

function counts_remove(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2})
    countsdelta = -get_counts(row, col, compsum, obsidxCounts)
    return countsdelta
end

function counts_rowswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2})
    countsdelta = get_counts(row, col, compsum, obsidxCounts) - get_counts(C.col2row[col], col, compsum, obsidxCounts)
    return countsdelta
end

function counts_colswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2})
    countsdelta = get_counts(row, col, compsum, obsidxCounts) - get_counts(row, C.row2col[row], compsum, obsidxCounts)
    return countsdelta
end

function counts_doubleswitch(row::Integer, col::Integer, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2})
    rowalt = C.col2row[col]
    colalt = C.row2col[row]
    countsdelta = get_counts(row, col, compsum, obsidxCounts)
    countsdelta += get_counts(rowalt, colalt, compsum, obsidxCounts)
    countsdelta -= get_counts(row, colalt, compsum, obsidxCounts)
    countsdelta -= get_counts(rowalt, col, compsum, obsidxCounts)
    return countsdelta
end

"""
determine type of move and returns change in number of links, block links, and stage links
"""
function pairmove_linkdelta(newrow::G, newcol::G, newblock::G, newstage::G, C::StageLinkMatrix{G}) where G <: Integer
    if iszero(C.row2col[newrow])
        if iszero(C.col2row[newcol]) #add link at newrow, newcol
            return one(G), newblock, newstage
        else #switch link from col2row[newcol], newcol -> newrow, newcol
            return zero(G), zero(G), zero(G)
        end
    elseif iszero(C.col2row[newcol]) #row not zero, col zero switch link from newrow, row2col[newrow] -> newrow, newcol
        return  zero(G), zero(G), zero(G)
    elseif C.row2col[newrow] == newcol #remove link from newrow, newcol
        return -one(G), -newblock, -newstage
    else
        return zero(G), zero(G), zero(G)
    end      
end

"""
    pairmove_loglikpCratio(newrow::G, newcol::G, C::StageLinkMatrix{G}, compsum::ComparisonSummary, logpCRatioAdd::T, logpCRatioRemove::T, loglikRatios::Array{T, 1}, loglikMissing::T = -Inf) where {G <: Integer, T <: AbstractFloat}

assume block and stage are known, calculate change in posterior, conditional on likelihood parameters, prior margins can be computed using lblockpriorratio
"""
function pairmove_loglikpCratio(row::G, col::G, C::StageLinkMatrix{G}, compsum::Union{ComparisonSummary, SparseComparisonSummary}, logpCRatioAdd::T, logpCRatioRemove::T, loglikRatios::Array{T, 1}, loglikMissing::T = -T(Inf)) where {G <: Integer, T <: AbstractFloat}

    if iszero(C.row2col[row])
        if iszero(C.col2row[col]) #add link at newrow, newcol
            return loglik_add(row, col, C, compsum, loglikRatios, loglikMissing) + logpCRatioAdd, false
        else #switch link from col2row[newcol], newcol -> newrow, newcol
            return loglik_rowswitch(row, col, C, compsum, loglikRatios, loglikMissing), false
        end
    elseif iszero(C.col2row[col]) #row not zero, col zero switch link from newrow, row2col[newrow] -> newrow, newcol
        return loglik_colswitch(row, col, C, compsum, loglikRatios, loglikMissing), false
    elseif C.row2col[row] == col #remove link from newrow, newcol
        return loglik_remove(row, col, C, compsum, loglikRatios, loglikMissing) + logpCRatioRemove, false
    else
        return loglik_doubleswitch(row, col, C, compsum, loglikRatios, loglikMissing), true
    end      
end

"""
    pairmove_countdelta(newrow::T, newcol::T, C::StageLinkMatrix{G}, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2}) where G <: Integer -> countdelta::Array{<: Integer, 1}

Sum changes in counts from pairmove
"""
function pairmove_countdelta(newrow::G, newcol::G, C::StageLinkMatrix{G}, compsum::Union{ComparisonSummary, SparseComparisonSummary}, obsidxCounts::Array{<:Integer, 2}) where G <: Integer
    countdelta = get_counts(newrow, newcol, compsum, obsidxCounts)
    if iszero(C.row2col[newrow])
        if iszero(C.col2row[newcol]) #add link at newrow, newcol
            return countdelta
        else #switch link from col2row[newcol], newcol -> newrow, newcol
            countdelta -= get_counts(C.col2row[newcol], newcol, compsum, obsidxCounts)
            return countdelta
        end
    elseif iszero(C.col2row[newcol]) #row not zero, col zero switch link from newrow, row2col[newrow] -> newrow, newcol
        countdelta -= get_counts(newrow, C.row2col[newrow], compsum, obsidxCounts)
        return countdelta
    elseif C.row2col[newrow] == newcol #remove link from newrow, newcol
        return -countdelta
    else
        countdelta  += get_counts(C.col2row[newcol], C.row2col[newrow], compsum, obsidxCounts)
        countdelta -= get_counts(C.col2row[newcol], newcol, compsum, obsidxCounts)
        countdelta -= get_counts(newrow, C.row2col[newrow], compsum, obsidxCounts)
        return countdelta
    end      
end

"""
    pairmove_update!(moveblock::G, movestage::G, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2}, loglikMissing::T = -Inf) where {G <: Integer, T <: AbstractFloat}

Randomwalk on record pairs contained in `moveblock` that are either unlinked or linked within the block.
"""
function pairmove_update!(moveblock::G, movestage::G, C::StageLinkMatrix,
                        compsum::Union{ComparisonSummary, SparseComparisonSummary},
                        bstages::BlockStages{G},
                        loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2},
                        logpCRatio::Function, loglikMissing::T = -Inf) where {G <: Integer, T <: AbstractFloat}

    rows = filter(row -> iszero(C.row2col[row]) || C.row2stage[row] == movestage, bstages.b2rows[moveblock])
    cols = filter(col -> iszero(C.col2row[col]) || C.row2stage[C.col2row[col]] == movestage, bstages.b2cols[moveblock])
    
    if length(rows) == 0 || length(cols) == 0
        return zeros(eltype(obsidxCounts), compsum.ncomp), false
    end

    ##Sample move
    moverow = sample(rows)
    movecol = sample(cols)

    #Resample if missing recordpair sampled
    if iszero(compsum.obsidx[row, col])
        return pairmove_update!(moveblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpCRatio, loglikMissing)
    end
    
    ##Calculate likelihoodratio
    logpCRatioAdd, logpCRatioRemove = logpCRatio(moveblock, C, bstages)
    lpratio, doubleswitch = pairmove_loglikpCratio(moverow, movecol, C, compsum, logpCRatioAdd, logpCRatioRemove, loglikRatios, loglikMissing)

    ##Re-sample p = 0.5 for switch moves since they are can sampled in two ways
    if doubleswitch
        if rand() < 0.5
            pairmove_update!(moveblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpCRatio, loglikMissing)
        end
    end
    
    countdelta = pairmove_countdelta(moverow, movecol, C, compsum, obsidxCounts)
    
    ##Perform move (or not)
    if rand() < exp(lpratio)
        C = pairmove_link!(moverow, movecol, movestage, C, bstages)
        return countdelta, true
    else
        return countdelta, false
    end
end

"""
log(sqrt(posterior(new) / posterior(old))) - assumes sqrt balancing function
"""
function log_move_weights(rows::Array{G, 1}, cols::Array{G, 1}, moveblock::G, C::StageLinkMatrix{G},
                          loglikRatios::Array{T, 1}, compsum::Union{ComparisonSummary, SparseComparisonSummary},
                          bstages::BlockStages{G},
                          logpCRatio::Function, log_balance_function::Function = lidentity_balance,
                          loglikMissing::T = -Inf) where {G <: Integer, T <: AbstractFloat}
    nr = length(rows)
    nc = length(cols)
    lmoveweights = fill(-Inf, nr * nc)
    logpCRatioAdd, logpCRatioRemove = logpCRatio(moveblock, C, bstages)
    idx = 0
    for cidx in 1:nc, ridx in 1:nr
        idx += 1
        lpratio, doubleswitch = pairmove_loglikpCratio(rows[ridx], cols[cidx], C, compsum, logpCRatioAdd, logpCRatioRemove, loglikRatios, loglikMissing)
        
        if doubleswitch
            lmoveweights[idx] = log_balance_function(lpratio) + loghalf
        else
            lmoveweights[idx] = log_balance_function(lpratio)
        end
    end
    lsummw = logsumexp(lmoveweights)
    return lmoveweights, lsummw
end

"""
    pairmove_locally_balanced_update!(moveblock::G, movestage::G, C::StageLinkMatrix, compsum::Union{ComparisonSummary, SparseComparisonSummary}, loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2}, loglikMissing::T = -Inf) where {G <: Integer, T <: AbstractFloat}

Locally balanced Randomwalk on record pairs contained in `moveblock` that are either unlinked or linked within the block.
"""
function pairmove_locally_balanced_update!(rows::Array{<:Integer, 1}, cols::Array{<:Integer, 1}, moveblock::G, movestage::G, C::StageLinkMatrix,
                                compsum::Union{ComparisonSummary, SparseComparisonSummary},
                                bstages::BlockStages{G},
                                loglikRatios::Array{T, 1}, obsidxCounts::Array{<:Integer, 2},
                                logpCRatio::Function, log_balance_function::Function = lidentity_balance,
                                loglikMissing::T = -T(Inf)) where {G <: Integer, T <: AbstractFloat}

    #rows = filter(row -> iszero(C.row2col[row]) || C.row2stage[row] == movestage, bstages.b2rows[moveblock])
    #cols = filter(col -> iszero(C.col2row[col]) || C.row2stage[C.col2row[col]] == movestage, bstages.b2cols[moveblock])

    if length(rows) == 0 || length(cols) == 0
        #println(moveblock)
        #println(movestage)
        #println(bstages)
        #println(C)
        #error("empty rows or columns")
        countdelta = zeros(eltype(obsidxCounts), size(obsidxCounts, 1))
        return C, countdelta, false
    end
    #println("Step")
    #println(moveblock)
    
    ##Move weights
    lmoveweights, lsummw = log_move_weights(rows, cols, moveblock, C, loglikRatios, compsum, bstages, logpCRatio, log_balance_function, loglikMissing)
    moverow, movecol, lmove = sample_proposal_full(rows, cols, lmoveweights, lsummw)

    #Can sample only missing entries if others are linked in different stages
    if iszero(compsum.obsidx[moverow, movecol])
        countdelta = zeros(eltype(obsidxCounts), size(obsidxCounts, 1))
        return C, countdelta, false
    end
    ##Sample move and reverse move
    #moveidx = sample(Weights(exp.(lmoveweights .- lsummw)))
    #println(moveidx)
    #movecidx, moveridx = divrem(moveidx, length(rows))

    #if moveridx == 0
    #    moveridx = length(rows)
    #else
    #    movecidx += 1
    #end
    
    #println(moveridx, ", ", movecidx)
    #moverow = rows[moveridx]
    #movecol = cols[movecidx]

    #if iszero(compsum.obsidx[moverow, movecol])
    #    @warn "Missing likelihood sampled, check code for bugs..."
    #end

    #println(moverow, ", ", movecol)
    #moverow, movecol = ind2sub((length(rows), length(cols)), moveidx)
    revrow, revcol = pairmove_inverse(moverow, movecol, movestage, C, bstages)
    #println(revrow, ", ", revcol)
    #if revrow == moverow
    #    revridx = moveridx
    #else
    #    revridx = findfirst(x -> x == revrow, rows)
    #end
    #if revcol == movecol
    #    revcidx = movecidx
    #else
    #    revcidx = findfirst(x -> x == revcol, cols)
    #end
    
    #revidx = (revcidx - 1) * length(rows) + revridx
    #println(revidx)
    #revidx = sub2ind((length(rows), length(cols)), revrow, revcol)
    
    #countdelta = pairmove_countdelta(moverow, movecol, C, compsum, obsidxCounts)
    logpCRatioAdd, logpCRatioRemove = logpCRatio(moveblock, C, bstages)
    loglikpCratio, doubleswitch = pairmove_loglikpCratio(moverow, movecol, C, compsum, logpCRatioAdd, logpCRatioRemove, loglikRatios, loglikMissing)
    if doubleswitch
        linvmove = log_balance_function(-loglikpCratio) + loghalf
    else
        linvmove = log_balance_function(-loglikpCratio)
    end
    
    ##Perform move
    C = pairmove_link!(moverow, movecol, movestage, C, bstages)
    
    ##Proposal weights and reverse move
    lpropweights, lsumpw = log_move_weights(rows, cols, moveblock, C, loglikRatios, compsum, bstages, logpCRatio, log_balance_function, loglikMissing)

    ##Compute ratio
    lmoveratio = linvmove - lmove + lsummw - lsumpw

    #also return change in likelihood
    if rand() < exp(loglikpCratio + lmoveratio)
        countdelta = -pairmove_countdelta(revrow, revcol, C, compsum, obsidxCounts)
        return C, countdelta, true
    else
        pairmove_link!(revrow, revcol, movestage, C, bstages)
        countdelta = zeros(eltype(obsidxCounts), size(obsidxCounts, 1))
        return C, countdelta, false
    end
end

"""
    stage_gibbs!(row::G, col::G, lpratio::T, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::Array{T, 1}, beta::Array{T, 1}) where {G <: Integer, T <: AbstractFloat}

Perform gibbs update on the link stage of a record pair.  `row` and `col` must either be linked to eachother or both be unlinked.
"""
function stage_gibbs!(row::G, col::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, logpCRatio::Function; minstage::Integer = one(G), maxstage::Integer = bstages.nstage) where {G <: Integer}
    if !iszero(C.row2col[row]) && C.row2col[row] != col
        @warn "Row link not to col, no update performed"
        return false, false
    elseif !iszero(C.col2row[col]) && C.col2row[col] != row
        @warn "Col link not to roq, no update performed"
        return false, false
    end

    if C.row2col[row] == col
        linkstart = true
        stagestart = C.row2stage[row]
        if iszero(row)
            @warn "row is zero??"
        end

        if iszero(col)
            @warn "col is zero??"
        end
        remove_link!(row, col, C, bstages)
    else
        linkstart = false
        stagestart = bstages.nstage + one(G)
    end

    logprop = zeros(Float64, bstages.nstage)
    for stage in one(G):bstages.nstage
        if bstages.rs2b[row, stage] == bstages.cs2b[col, stage] && !iszero(bstages.rs2b[row, stage]) && stage >= minstage && stage <= maxstage
            logprop[stage] = logpCRatio(bstages.rs2b[row, stage], C, bstages)[1]
        else
            logprop[stage] = -Inf
        end
    end
    #logprop[end] = -lpratio
    softmax!(logprop)
    updatestage = G(sample(Weights(logprop)))
    if updatestage <= bstages.nstage
        add_link!(row, col, updatestage, C, bstages)
        linkmove = !linkstart
    else
        linkmove = linkstart
    end
    
    return linkmove, stagestart != updatestage
end

function pairmove_globally_balanced_update!(rows::Array{<:Integer, 1}, cols::Array{<:Integer, 1}, moveblock::G, movestage::G, C::StageLinkMatrix,
                                 compsum::Union{ComparisonSummary, SparseComparisonSummary},
                                 bstages::BlockStages{G},
                                 loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2},
                                 logpCRatio::Function) where G <: Integer
    return pairmove_locally_balanced_update!(rows, cols, moveblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpCRatio, identity)
end

function pairmove_locally_balanced_sqrt_update!(rows::Array{<:Integer, 1}, cols::Array{<:Integer, 1}, moveblock::G, movestage::G, C::StageLinkMatrix,
                                     compsum::Union{ComparisonSummary, SparseComparisonSummary},
                                     bstages::BlockStages{G},
                                     loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2},
                                     logpCRatio::Function) where G <: Integer
    return pairmove_locally_balanced_update!(rows, cols, moveblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpCRatio, lsqrt)
end

function pairmove_locally_balanced_barker_update!(rows::Array{<:Integer, 1}, cols::Array{<:Integer, 1}, moveblock::G, movestage::G, C::StageLinkMatrix,
                                       compsum::Union{ComparisonSummary, SparseComparisonSummary},
                                       bstages::BlockStages{G},
                                       loglikRatios::Array{<:AbstractFloat, 1}, obsidxCounts::Array{<:Integer, 2},
                                       logpCRatio::Function) where G <: Integer
    return pairmove_locally_balanced_update!(rows, cols, moveblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpCRatio, lbarker)
end

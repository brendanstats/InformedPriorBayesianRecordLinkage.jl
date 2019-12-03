#using SpecialFunctions
#using StatsFuns
#using StatsBase
#using Random

#include("blockstages.jl")
#include("stagelinkmatrix.jl")
#w = weights()
#sample(1:nstage, w)
#lfactorial
"""
Evaluate block beta bipartite prior (one stage)
"""
function blockbetabipartite(blocklinks::Array{<:Integer, 1}, blockmindim::Array{<:Integer, 1},
                            alpha::T, beta::T,
                            totblocklinks::Integer = sum(blocklinks),
                            totmindim::Integer = sum(blockmindim)) where {T <: AbstractFloat}
    c = zero(T)
    for (l, r) in zip(blocklinks, blockmindim)
        if r > l
            c ./ prod((r - l + 1):l)
        end
    end
    return c * beta(totblocklinks + alpha, totmindim - totblocklinks + beta) / beta(alpha, beta)
end

"""
Evaluate log(block beta bipartite prior) (one stage)
"""
function lblockbetabipartite(blocklinks::Array{<:Integer, 1}, blockmindim::Array{<:Integer, 1},
                             alpha::T, beta::T,
                             totblocklinks::Integer = sum(blocklinks),
                             totmindim::Integer = sum(blockmindim)) where {T <: AbstractFloat}
    lc = zero(T)
    for (l, r) in zip(blocklinks, blockmindim)
        lc += (lfactorial(r - l) - lfactorial(r))
    end
    return lc + lbeta(totblocklinks + alpha, totmindim - totblocklinks + beta) - lbeta(alpha, beta)
end

"""
Evaluate ratio of block beta bipartite prior(links1) / block beta bipartite prior(links2) (one stage)
"""
function blockbetabipartite_ratio(blocklinks1::Array{<:Integer, 1}, blocklinks2::Array{<:Integer, 1},
                                  blockmindim::Array{<:Integer, 1}, alpha::T, beta::T,
                                  totblocklinks1::Integer = sum(blocklinks1),
                                  totmindim::Integer = sum(blockmindim)) where {T <: AbstractFloat}
    c = zero(T)
    for (l1, l2, r) in zip(blocklinks1, blocklinks2, blockmindim)
        c *= factorial(blocklinks1 - blockmindim) / factorial(blocklinks2 - blockmindim)
    end
    return c * beta(totblocklinks1 + alpha, totmindim - totblocklinks1 + beta) / beta(totblocklinks2 + alpha, totmindim - totblocklinks2 + beta)
end

"""
Evaluate log ratio of block beta bipartite prior(links1) / block beta bipartite prior(links2) (one stage)
"""
function lblockbetabipartite_ratio(blocklinks1::Array{<:Integer, 1}, blocklinks2::Array{<:Integer, 1},
                                   blockmindim::Array{<:Integer, 1}, alpha::T, beta::T,
                                   totblocklinks1::Integer = sum(blocklinks1),
                                   totmindim::Integer = sum(blockmindim)) where {T <: AbstractFloat}
    lc = zero(T)
    for (l1, l2, r) in zip(blocklinks1, blocklinks2, blockmindim)
        if l1 > l2
            lc += lfactorial(l1 - r) - lfactorial(l2 - r)
        elseif l1 < l2
            lc += lfactorial(l2 - r) - lfactorial(l1 - r)
        end
    end
    return lc + lbeta(totblocklinks1 + alpha, totmindim - totblocklinks1 + beta) - lbeta(totblocklinks2 + alpha, totmindim - totblocklinks2 + beta)
end

"""
Evaluate log ratio of block beta bipartite prior(links1) / block beta bipartite prior(links2) (multiple stage)
"""
function lblockbetabipartite_sequence(blocklinks::Array{<:Integer, 1},
                                      blockprevlinks::Array{<:Integer, 1},
                                      blockmindim::Array{<:Integer, 1},
                                      stagelinks::Array{<:Integer, 1},
                                      stageprevlinks::Array{<:Integer, 1},
                                      stagemindim::Array{<:Integer, 1},
                                      stagealpha::Array{T, 1}, stagebeta::Array{T, 1},
                                      nblocks::Integer = length(blocklinks),
                                      nstage::Integer = length(stagealpha)) where {T <: AbstractFloat}
    lprior = zero(T)
    for (lb, pl, bmin) in zip(blocklinks, blockprevlinks, blockmindim)
        r = bmin - pl
        lprior += lfactorial(r - lb) - lfactorial(r)
    end

    for (sl, cpsl, smin, alpha, beta) in zip(stagelinks, stageprevlinks, stagemindim, stagealpha, stagebeta)
        lprior += lbeta(sl + alpha, smin - cpsl - sl + beta) - lbeta(alpha, beta)
    end
    return lprior
end

#(log(prior(add link to stage) / prior(current)
function lblockbetabipartite_sequencep1(block::Integer,
                                        stage::Integer,
                                        blocklinks::Array{<:Integer, 1},
                                        blockprevlinks::Array{<:Integer, 1},
                                        blockmindim::Array{<:Integer, 1},
                                        stagelinks::Array{<:Integer, 1},
                                        stageprevlinks::Array{<:Integer, 1},
                                        stagemindim::Array{<:Integer, 1},
                                        stagealpha::Array{T, 1}, stagebeta::Array{T, 1},
                                        nblocks::Integer = length(blocklinks),
                                        nstage::Integer = length(stagealpha)) where {T <: AbstractFloat}
    lprior = zero(T)
    for (lb1, lb2, pl1, pl2, bmin) in zip(blocklinks1, blocklinks2, blockprevlinks1, blockprevlinks2, blockmindim)
        r1 = bmin - pl1
        r2 = bmin - pl2
        lprior += lfactorial(r1 - lb1) - lfactorial(r1) - lfactorial(r1 - lb1) + lfactorial(r1)
    end

    for (sl, cpsl, smin, alpha, beta) in zip(stagelinks, stageprevlinks, stagemindim, stagealpha, stagebeta)
        lprior += lbeta(sl + alpha, smin - cpsl - sl + beta) - lbeta(alpha, beta)
    end
    return lprior
end

"""
    lp1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T) -> logpriorratio
    lp1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat} -> logpriorratio

Return log block beta-bipartite ratio for adding link to a block: log(add/current)
"""
function lp1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    if blocklinks > condblockmin || blocklinks < zero(G)
        return -Inf
    end
    return log(blocklinks + alpha) - log(condstagemin - blocklinks - one(G) + beta) - log(condblockmax - blocklinks)
end

function lp1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    stage = bstages.b2s[block]
    return lp1(C.nlinkblock[block], C.nlinkstage[stage], C.block2condmin[block], C.block2condmax[block], C.stage2condmin[stage], alpha, beta)
end

"""
    lm1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T) -> logpriorratio
    lm1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat} -> logpriorratio

Return log block beta-bipartite ratio for removing link from block: log(remove/current)
"""
function lm1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    if blocklinks > condblockmin || blocklinks < zero(G)
        return -Inf
    end
    return log(condstagemin - blocklinks + one(G) + beta) - log(blocklinks - one(G) + alpha) + log(condblockmax - blocklinks + one(G))
end

function lm1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    stage = bstages.b2s[block]
    return lm1(C.nlinkblock[block], C.nlinkstage[stage], C.block2condmin[block], C.block2condmax[block], C.stage2condmin[stage], alpha, beta)
end

"""
    le1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T) -> logpriorratio
    le1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat} -> logpriorratio

Return log block beta-bipartite ratio for expanding block by one: log(expanded/current)
"""
function le1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    if blocklinks > condblockmin || blocklinks < zero(G)
        return -Inf
    end
    return log(condblockmax - blocklinks + one(G)) - log(condblockmax + one(G)) + log(condstagemin - blocklinks + one(G) + beta) - log(condstagemin + alpha + beta)
end

function le1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    stage = bstages.b2s[block]
    return le1(C.nlinkblock[block], C.nlinkstage[stage], C.block2condmin[block], C.block2condmax[block], C.stage2condmin[stage], alpha, beta)
end

"""
    ls1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T) -> logpriorratio
    ls1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat} -> logpriorratio

Return log block beta-bipartite ratio for shrinking block by one: log(shrunk/current)
"""
function ls1(blocklinks::G, stagelinks::G, condblockmin::G, condblockmax::G, condstagemin::G, alpha::T, beta::T)  where {G <: Integer, T <: AbstractFloat}
    if blocklinks > condblockmin || blocklinks < zero(G)
        return -Inf
    end
    return log(condblockmax) - log(condblockmax - blocklinks) + log(condstagemin - one(G) + alpha + beta) - log(condstagemin - one(G) - blocklinks + beta)
end

function ls1(block::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T) where {G <: Integer, T <: AbstractFloat}
    stage = bstages.b2s[block]
    return ls1(C.nlinkblock[block], C.nlinkstage[stage], C.block2condmin[block], C.block2condmax[block], C.stage2condmin[stage], alpha, beta)
end

"""
   lblockpriorratio(moveblock::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::T, beta::T) -> addlogprior, removelogprior

Compute the stage block beta-bipartite prior ratio for both adding and remove link to `moveblock`, returns log(add / current), log(remove, current)
"""
function lblockpriorratio(moveblock::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}, alpha::Array{T, 1}, beta::Array{T, 1}) where {G <: Integer, T <: AbstractFloat}
    addlp = lp1(moveblock, C, bstages, alpha[bstages.b2s[moveblock]], beta[bstages.b2s[moveblock]])
    removelp = lm1(moveblock, C, bstages, alpha[bstages.b2s[moveblock]], beta[bstages.b2s[moveblock]])
    for block in bstages.b2conditionedb[moveblock]
        addlp += ls1(block, C, bstages, alpha[bstages.b2s[block]], beta[bstages.b2s[block]])
        removelp += le1(block, C, bstages, alpha[bstages.b2s[block]], beta[bstages.b2s[block]])
    end
    return addlp, removelp
end

function draw_stageblockbetabipartite(bstages::BlockStages, alphas::Array{<:AbstractFloat, 1}, betas::Array{<:AbstractFloat, 1})
    if length(alphas) != length(betas)
        error("alphas and betas must be the same length")
    elseif length(alphas) != bstages.nstage
        error("parameters must be the same length as the number of stages")
    end
    
    C = StageLinkMatrix(bstages)
    for stage in 1:bstages.nstage
        p = rand(Beta(alphas[stage], betas[stage]))
        for block in findall(bstages.b2s .== stage)
            if C.block2condmin[block] > 0
                n = rand(Binomial(C.block2condmin[block], p))
                rows = sample(get_blockopenrows(block, C, bstages), n, replace = false)
                cols = sample(get_blockopencols(block, C, bstages), n, replace = false)
                for (row, col) in zip(rows, cols)
                    add_link!(row, col, stage, C, bstages)
                end
            end
        end
    end
    return C
end

#permutedims(mapreduce(x -> firststage_counts(draw_stageblockbetabipartite(bstages, [1.0,1.0], [1.0,1.0]), bstages), hcat, 1:1000))

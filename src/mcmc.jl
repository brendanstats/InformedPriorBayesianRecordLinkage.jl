function counts_matches(C::StageLinkMatrix{G}, compsum::Union{ComparisonSummary, SparseComparisonSummary}) where G <: Integer

    #count occurences of each observation in obsvecs
    matchvecct = zeros(Int64, length(compsum.obsvecct))
    if compsum.nrow < compsum.ncol
        for (ii, jj) in pairs(IndexLinear(), C.row2col)
            if jj != zero(G)
                matchvecct[compsum.obsidx[ii, jj]] += 1
            end
        end
    else
        for (jj, ii) in pairs(IndexLinear(), C.col2row)
            if ii != zero(G)
                matchvecct[compsum.obsidx[ii, jj]] += 1
            end
        end
    end

    #map observation occurences to counts
    matchcounts = zeros(Int64, length(compsum.counts))
    matchobs = zeros(Int64, compsum.ncomp)

    for (jj, ct) in pairs(IndexLinear(), matchvecct)
        if ct > 0
            for ii in 1:compsum.ncomp
                if compsum.obsvecs[ii, jj] != 0
                    matchobs[ii] += ct
                    matchcounts[compsum.cadj[ii] + compsum.obsvecs[ii, jj]] += ct
                end
            end
        end
    end
    return matchcounts, matchobs
end

"""
bstages = BlockStages([[1, 1, 2, 2] [3, 3, 3, 3]], [[1, 1, 2, 2] [3, 3, 3, 3]])
phb = PosthocBlocks(Dict(1 => [1], 2 => [2, 3], 3 => [4]), Dict(1 => [1], 2 => [2, 3], 3 => [4]), [1, 2, 1], [1, 2, 1], [true, false, true], [1, 4, 1], 4, 4, 3, 6)
restrict_blocks(phb, bstages)
Dicts mapping (posthoc block, stage block) -> rows, (posthoc block, stage block) -> cols 
"""
function restrict_blocks(phb::PosthocBlocks{G}, bstages::BlockStages{T}) where {G <: Integer, T <: Integer}
    
    rowcounts = spzeros(bstages.nblock, phb.nblock)
    colcounts = spzeros(bstages.nblock, phb.nblock)
    for pblock in one(G):phb.nblock
        rows = phb.block2rows[pblock]
        cols = phb.block2cols[pblock]
        for stage in one(T):bstages.nstage
            for row in rows
                if !iszero(bstages.rs2b[row, stage])
                    rowcounts[bstages.rs2b[row, stage], pblock] += 1
                end
            end

            for col in cols
                if !iszero(bstages.cs2b[col, stage])
                    colcounts[bstages.cs2b[col, stage], pblock] += 1
                end
            end
        end
    end

    restrictedRows = Dict{Tuple{G, T}, Array{G, 1}}()
    restrictedCols = Dict{Tuple{G, T}, Array{G, 1}}()
    
    sblocks = rowvals(rowcounts)
    for pblock in one(G):phb.nblock
        rows = phb.block2rows[pblock]
        cols = phb.block2cols[pblock]
        for idx in nzrange(rowcounts, pblock)
            stageblock = sblocks[idx]
            if colcounts[stageblock, pblock] > 0
                stage = bstages.b2s[stageblock]
                restrictedRows[(pblock, T(stageblock))] = rows[bstages.rs2b[rows, stage] .== stageblock]
                restrictedCols[(pblock, T(stageblock))] = cols[bstages.cs2b[cols, stage] .== stageblock]
            end
        end
    end

    return restrictedRows, restrictedCols
end

function filter_rows_stage(rows::Array{<:Integer, 1}, C::StageLinkMatrix, stage::Integer)
    return filter(row -> iszero(C.row2col[row]) || C.row2stage[row] == stage, rows)
end

function filter_cols_stage(cols::Array{<:Integer, 1}, C::StageLinkMatrix, stage::Integer)
    return filter(col -> iszero(C.col2row[col]) || C.row2stage[C.col2row[col]] == stage, cols)
end

function mh_gibbs_count(
    nsteps::Integer,
    C0::StageLinkMatrix{G},
    compsum::Union{ComparisonSummary, SparseComparisonSummary},
    phb::PosthocBlocks,
    bstages::BlockStages{G},
    priorM::Array{T, 1},
    priorU::Array{T, 1},
    logpdfC::Function,
    transitionLinkC!::Function,
    transitionStageC!::Function;
    minstage::Integer = one(G),
    maxstage::Integer = bstages.nstage,
    mingibbsstage::Integer = one(G),
    maxgibbsstage::Integer = bstages.nstage) where {G <: Integer, T <: Real}

    restrictedRows, restrictedCols = restrict_blocks(phb, bstages)
    
    #MCMC Chains
    CArray = zeros(Int64, bstages.nrow, bstages.ncol, bstages.nstage)
    nlinkArray = zeros(G, bstages.nstage, nsteps)
    MArray = Array{Float64}(undef, length(priorM), nsteps)
    UArray = Array{Float64}(undef, length(priorU), nsteps)
    transLink = zero(Int64)
    transStage = zero(Int64)
    
    ##Initial States
    obsidxCounts = get_obsidxcounts(compsum) #each column is an observation
    #obsDeltas = obs_delta(compsum)
    C = deepcopy(C0)
    matchcounts, matchobs = counts_matches(C, compsum)
    pM, pU, loglikRatios = gibbs_MU_draw(matchcounts, compsum, obsidxCounts, priorM, priorU)
    
    #Outer iteration (recorded)
    for ii in 1:nsteps
        
            for (pblock, sblock) in keys(restrictedRows)
                #only attempt move if block can be sampled
                if C.block2condmin[sblock] > zero(G) && bstages.b2s[sblock] >= minstage && bstages.b2s[sblock] <= maxstage
                    movestage = bstages.b2s[sblock]
                    rows = filter_rows_stage(restrictedRows[(pblock, sblock)], C, movestage)
                    cols = filter_cols_stage(restrictedCols[(pblock, sblock)], C, movestage)
                    #println(rows)
                    C, countdelta, move = transitionLinkC!(rows, cols, sblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpdfC)
                    
                    if move
                        transLink += one(Int64)
                        matchcounts += countdelta
                    end
                end
            end

        ##Update link stages
        for row in one(G):bstages.nrow
            if !iszero(C.row2col[row]) && minstage < maxstage && C.row2stage[row] >= minstage && C.row2stage[row] <= maxstage
                countdelta = pairmove_countdelta(row, C.row2col[row], C, compsum, obsidxCounts)
                linkmove, stagemove = transitionStageC!(row, C.row2col[row], C, bstages, logpdfC, minstage = mingibbsstage, maxstage = maxgibbsstage)
                if linkmove
                    matchcounts += countdelta
                    transLink += one(Int64)
                end
                if stagemove
                    transStage += one(Int64)
                end
            end
        end
        
        ##Update parameters
        pM, pU, loglikRatios = gibbs_MU_draw(matchcounts, compsum, obsidxCounts, priorM, priorU)
        
        #Add states to chain
        for row in one(G):bstages.nrow
            if !iszero(C.row2col[row])
                CArray[row, C.row2col[row], C.row2stage[row]] += one(Int64)
            end
        end
        nlinkArray[:, ii] = C.nlinkstage   
        MArray[:, ii] = pM
        UArray[:, ii] = pU
    end
    
    return ParameterChain(counts2indicies(CArray), permutedims(nlinkArray), permutedims(MArray), permutedims(UArray), nsteps, false), transLink, transStage, C

end

function mh_gibbs_trace(
    nsteps::Integer,
    C0::StageLinkMatrix{G},
    compsum::Union{ComparisonSummary, SparseComparisonSummary},
    phb::PosthocBlocks,
    bstages::BlockStages{G},
    priorM::Array{T, 1},
    priorU::Array{T, 1},
    logpdfC::Function,
    transitionLinkC!::Function,
    transitionStageC!::Function;
    minstage::Integer = one(G),
    maxstage::Integer = bstages.nstage,
    mingibbsstage::Integer = one(G),
    maxgibbsstage::Integer = bstages.nstage) where {G <: Integer, T <: Real}

    restrictedRows, restrictedCols = restrict_blocks(phb, bstages)
    
    #MCMC Chains
    outrows = Int[]
    outcols = Int[]
    outstages = Int[]
    outstart = Int[]
    outstop = Int[]
    nlinkArray = zeros(G, bstages.nstage, nsteps)
    MArray = Array{Float64}(undef, length(priorM), nsteps)
    UArray = Array{Float64}(undef, length(priorU), nsteps)
    transLink = zero(Int64)
    transStage = zero(Int64)

    ##Initial States
    obsidxCounts = get_obsidxcounts(compsum)
    C = deepcopy(C0)
    matchcounts, matchobs = counts_matches(C, compsum)
    pM, pU, loglikRatios = gibbs_MU_draw(matchcounts, compsum, obsidxCounts, priorM, priorU)

    currrow2col = copy(C.row2col)
    currrow2stage = copy(C.row2stage)
    startrow2col = ones(Int, length(currrow2col)) #zeros(Int, length(currrow2col))
    
    #Outer iteration (recorded)
    for ii in 1:nsteps
        
        for (pblock, sblock) in keys(restrictedRows)
            #only attempt move if block can be sampled
            if C.block2condmin[sblock] > zero(G) && bstages.b2s[sblock] >= minstage && bstages.b2s[sblock] <= maxstage
                movestage = bstages.b2s[sblock]
                rows = filter_rows_stage(restrictedRows[(pblock, sblock)], C, movestage)
                cols = filter_cols_stage(restrictedCols[(pblock, sblock)], C, movestage)
                #println(rows)
                C, countdelta, move = transitionLinkC!(rows, cols, sblock, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpdfC)
                
                if move
                    transLink += one(Int64)
                    matchcounts += countdelta
                end
            end
        end
            
        ##Update link stages
        for row in one(G):bstages.nrow
            if !iszero(C.row2col[row]) && minstage < maxstage && C.row2stage[row] >= minstage && C.row2stage[row] <= maxstage
                countdelta = pairmove_countdelta(row, C.row2col[row], C, compsum, obsidxCounts)
                linkmove, stagemove = transitionStageC!(row, C.row2col[row], C, bstages, logpdfC, minstage = mingibbsstage, maxstage = maxgibbsstage)
                if linkmove
                    matchcounts += countdelta
                    transLink += one(Int64)
                end
                if stagemove
                    transStage += one(Int64)
                end
            end
        end
        
        ##Update parameters
        pM, pU, loglikRatios = gibbs_MU_draw(matchcounts, compsum, obsidxCounts, priorM, priorU)
        
        #Add states to chain
        for row in one(G):bstages.nrow
            if currrow2col[row] != C.row2col[row] || currrow2stage[row] != C.row2stage[row]
                #record if deletion or move
                if !iszero(currrow2col[row])
                    push!(outrows, row)
                    push!(outcols, currrow2col[row])
                    push!(outstages, currrow2stage[row])
                    push!(outstart, startrow2col[row])
                    push!(outstop, ii - 1)                    
                end
                
                currrow2col[row] = C.row2col[row]
                currrow2stage[row] = C.row2stage[row]
                startrow2col[row] = ii
            end
        end
        nlinkArray[:, ii] = C.nlinkstage
        MArray[:, ii] = pM
        UArray[:, ii] = pU
    end
    
    for row in one(G):bstages.nrow
        if !iszero(C.row2col[row])
            push!(outrows, row)
            push!(outcols, currrow2col[row])
            push!(outstages, currrow2stage[row])
            push!(outstart, startrow2col[row])
            push!(outstop, nsteps)                    
        end
    end
    return ParameterChain([outrows outcols outstages outstart outstop][outstart .<= outstop, :], permutedims(nlinkArray), permutedims(MArray), permutedims(UArray), nsteps, true), transLink, transStage, C
end

function mh_gibbs_count_fixedparam(
    nsteps::Integer,
    C0::StageLinkMatrix{G},
    compsum::Union{ComparisonSummary, SparseComparisonSummary},
    bstages::BlockStages{G},
    pM::Array{T, 1},
    pU::Array{T, 1},
    logpdfC::Function,
    transitionLinkC!::Function,
    transitionStageC!::Function;
    minstage::Integer = one(G),
    maxstage::Integer = bstages.nstage,
    mingibbsstage::Integer = one(G),
    maxgibbsstage::Integer = bstages.nstage) where {G <: Integer, T <: Real}

    #MCMC Chains
    CArray = zeros(Int64, bstages.nrow, bstages.ncol, bstages.nstage)
    nlinkArray = zeros(G, bstages.nstage, nsteps)
    transLink = zero(Int64)
    transStage = zero(Int64)

    ##Initial States
    obsidxCounts = get_obsidxcounts(compsum) #each column is an observation
    C = deepcopy(C0)
    logDiff = log.(pM) - log.(pU)
    loglikRatios = obsidxCounts' * logDiff

    #Outer iteration (recorded)
    for ii in 1:nsteps

        for bb in one(G):bstages.nblock
            if C.block2condmin[bb] > zero(G) && bstages.b2s[bb] >= minstage && bstages.b2s[bb] <= maxstage
                
                movestage = bstages.b2s[bb]
                C, countdelta, move = transitionLinkC!(bb, movestage, C, compsum, bstages, loglikRatios, obsidxCounts, logpdfC)
                
                if move
                    transLink += one(Int64)
                end
            end
        end
        
        ##Update link stages
        for row in one(G):bstages.nrow
            if !iszero(C.row2col[row]) && minstage < maxstage && C.row2stage[row] >= minstage && C.row2stage[row] <= maxstage
                
                linkmove, stagemove = transitionStageC!(row, C.row2col[row], C, bstages, logpdfC, minstage = mingibbsstage, maxstage = maxgibbsstage)
                
                if linkmove
                    transLink += one(Int64)
                end
                if stagemove
                    transStage += one(Int64)
                end
            end
        end
    end
    #Add states to chain
    for row in one(G):bstages.nrow
        if !iszero(C.row2col[row])
            CArray[row, C.row2col[row], C.row2stage[row]] += one(Int64)
        end
    end
    nlinkArray[:, ii] = C.nlinkstage        
    
    return CArray, permutedims(nlinkArray), transLink, transStage, C
end

#Add sequential function to fix at different thresholds
function mh_gibbs_trace_threshold(
    burnin::Integer,
    nsteps::Array{<:Integer, 1},
    threshold::Array{<:AbstractFloat, 1},
    C0::StageLinkMatrix{G},
    compsum::Union{ComparisonSummary, SparseComparisonSummary},
    bstages::BlockStages{G},
    priorM::Array{T, 1},
    priorU::Array{T, 1},
    logpdfC::Function,
    transitionLinkC!::Function,
    transitionStageC!::Function;
    startstage::Integer = one(G),
    endstage::Integer = bstages.nstage) where {G <: Integer, T <: Real}

    outChain = Dict{Int, ParameterChain{Int, T}}()
    outC = Dict{Int, typeof(C0)}()

    #Run burnin
    C = deepcopy(C0)
    pchain, linkmoves, stagemoves, C = mh_gibbs_trace_inplace(burnin, C, compsum, bstages, priorM, priorU,  logpdfC, transitionLinkC!, transitionStageC!, minstage = startstage, maxstage = startstage, mingibbsstage = startstage, maxgibbsstage = startstage)

    #Collect Samples
    for stage in startstage:endstage
        pchain, linkmoves, stagemoves, C = mh_gibbs_trace_inplace(nsteps[stage], C, compsum, bstages, priorM, priorU,  logpdfC, transitionLinkC!, transitionStageC!, minstage = stage, maxstage = stage, mingibbsstage = stage, maxgibbsstage = stage)
        outChain[stage] = pchain
        stagects = get_linkstagecounts(pchain)
        stagethreshold = floor(Int, threshold[stage] * nsteps[stage])
        if stage == startstage
            C = StageLinkMatrix(bstages)
        else
            C = deepcopy(outC[stage - 1])
        end
        for ii in 1:size(stagects, 1)
            if stagects[ii, end] > stagethreshold
                add_link!(stagects[ii, 1], stagects[ii, 2], stagects[ii, 3], C, bstages)
            end
        end
        outC[stage] = deepcopy(C)
    end
    
    return outChain, outC
end

#Add sequential function to take multiple samples
function mh_gibbs_trace_sequential_sample(
    burnin::Integer,
    nsteps::Array{<:Integer, 1},
    nsample::Array{<:Integer, 1},
    C0::StageLinkMatrix{G},
    compsum::Union{ComparisonSummary, SparseComparisonSummary},
    bstages::BlockStages{G},
    priorM::Array{T, 1},
    priorU::Array{T, 1},
    logpdfC::Function,
    transitionLinkC!::Function,
    transitionStageC!::Function;
    startstage::Integer = one(G),
    endstage::Integer = bstages.nstage) where {G <: Integer, T <: Real}

    pastC = Dict{Array{Int, 1}, typeof(C0)}()
    
    #Run burnin
    C = deepcopy(C0)
    pchain, linkmoves, stagemoves, C = mh_gibbs_trace_inplace(burnin, C, compsum, bstages, priorM, priorU,  logpdfC, transitionLinkC!, transitionStageC!, minstage = startstage, maxstage = startstage, mingibbsstage = startstage, maxgibbsstage = startstage)
    pastC[Int[]] = C
    #Collect Samples
    for stage in startstage:endstage
        nextC = Dict{Array{Int, 1}, typeof(C0)}()
        for (snum, sC) in pastC
            for kk in 1:nsample[stage]
                pchain, linkmoves, stagemoves, C = mh_gibbs_trace_inplace(nsteps[stage], sC, compsum, bstages, priorM, priorU,  logpdfC, transitionLinkC!, transitionStageC!, minstage = stage, maxstage = stage,  mingibbsstage = stage, maxgibbsstage = stage)
                nextC[[snum; kk]] = C
            end
        end
        pastC = nextC
    end
    return pastC
end

struct StageLinkMatrix{G <: Integer}
    row2col::Array{G, 1}
    row2block::Array{G, 1}
    row2stage::Array{G, 1}
    col2row::Array{G, 1}    
    nlinkblock::Array{G, 1}
    nlinkstage::Array{G, 1}
    block2condmin::Array{G, 1}
    block2condmax::Array{G, 1}
    stage2condmin::Array{G, 1}    
end

function StageLinkMatrix(nrow::G, ncol::G, bmindim::Array{G, 1}, bmaxdim::Array{G, 1}, smindim::Array{G, 1}) where {G <: Integer}
    return StageLinkMatrix(zeros(G, nrow), zeros(G, nrow), zeros(G, nrow), zeros(G, ncol), zeros(G, length(bmindim)), zeros(G, length(smindim)), copy(bmindim), copy(bmaxdim), copy(smindim))
end

StageLinkMatrix(bstages::BlockStages{G}) where {G <: Integer} = StageLinkMatrix(bstages.nrow, bstages.ncol, bstages.bmindim, bstages.bmaxdim, bstages.smindim)

function StageLinkMatrix(row2col::Array{G, 1}, row2stage::Array{G, 1}, bstages::BlockStages{G}) where {G <: Integer}
    if length(row2col) != bstages.nrow
        error("length of row2col and bstages.nrow must match")
    elseif length(row2col) != length(row2stage)
        error("length of row2col and row2stage  must match")
    end
    C = StageLinkMatrix(bstages)
    for row in one(G):bstages.nrow
        if !iszero(row2col[row])
            C = add_link!(row, row2col[row], row2stage[row], C, bstages)
        end
    end
    return C
end

"""
    add_link!(row::G, col::G, stage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer} -> C

Link `row` and `col` in `stage` checking that neither is already linked and stage is feasible and update conditonal block sizes.

See also: [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function add_link!(row::G, col::G, stage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer}
    if !iszero(C.row2col[row])
        @warn "row non-empty no addition made"
    elseif !iszero(C.col2row[col])
        @warn "col non-empty no addition made"
    elseif iszero(bstages.rs2b[row, stage])
        @warn "row, col, stage set not included in any block, no addition made"
    else
        block = bstages.rs2b[row, stage]
        C.row2col[row] = col
        C.row2block[row] = block
        C.row2stage[row] = stage
        C.col2row[col] = row
        for bb in bstages.b2conditionedb[C.row2block[row]]
            C.block2condmin[bb] -= one(G)
            C.block2condmax[bb] -= one(G)
        end
        for ss in (C.row2stage[row] + one(G)):bstages.nstage
            C.stage2condmin[ss] -= one(G)
        end
        C.nlinkblock[block] += one(G)
        C.nlinkstage[stage] += one(G)
    end
    return C
end

"""
    remove_link!(row::G, col::G, stage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer} -> C

Delete `row` and `col` link checking that link exists and update conditonal block sizes.

See also: [`add_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function remove_link!(row::G, col::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where {G <: Integer}
    if C.row2col[row] != col
        warn("Row and col not linked, no deletion made")
    else
        if iszero(C.row2block[row])
            println(C)
        end
        for block in bstages.b2conditionedb[C.row2block[row]]
            C.block2condmin[block] += one(G)
            C.block2condmax[block] += one(G)
        end
        for stage in (C.row2stage[row] + one(G)):bstages.nstage
            C.stage2condmin[stage] += one(G)
        end
        C.nlinkblock[C.row2block[row]] -= one(G)
        C.nlinkstage[C.row2stage[row]] -= one(G)
        C.row2col[row] = zero(G)
        C.row2block[row] = zero(G)
        C.row2stage[row] = zero(G)
        C.col2row[col] = zero(G)
    end
    return C
end

"""
    stageswitch_link!(row::G, col::G, newblock::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> C

Change stage of `row` and `col` link checking that link exists and update conditional block sizes.

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function stageswitch_link!(row::G, col::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(row)
        @warn "row must be non-zero"
    elseif iszero(col)
        @warn "col must be non-zero"
    elseif iszero(newstage)
        @warn "newstage must be non-zero"
    elseif C.row2col[row] != col
        @warn "Row and col not linked, no switch made"
    elseif bstages.rs2b[row, newstage] != bstages.cs2b[col, newstage]
        @warn "row and column not in same block in stage $newstage, no switch made"
    elseif iszero(bstages.rs2b[row, newstage])
        @warn "row not included in newstage $newstage, no switch made"
    elseif iszero(bstages.cs2b[col, newstage])
        @warn "col not included in newstage $newstage, no switch made"
    elseif bstages.rs2b[row, newstage] != bstages.cs2b[col, newstage]
        @warn "row and column not in the same stage in stage $newstage, no switch made"
    elseif  C.row2stage[row] != newstage
        block = C.row2block[row]
        stage = C.row2stage[row]
        newblock = bstages.rs2b[row, newstage]
        
        C.row2block[row] = newblock
        C.row2stage[row] = newstage
        C.nlinkblock[C.nlinkstage[row]] -= one(G)
        C.nlinkstage[C.row2stage[row]] -= one(G)
        C.nlinkblock[newblock] += one(G)
        C.nlinkstage[newstage] += one(G)

        if stage < newstage
            for bb in bstages.rs2b[(stage + one(G)):newstage]
                C.block2condmin[bb] += one(G)
                C.block2condmax[bb] += one(G)
            end
            for ss in (stage + one(G)):newstage
                C.stage2condmin[ss] += one(G)
            end
        elseif newstage < stage
            for bb in bstages.rs2b[(newstage + one(G)):stage]
                C.block2condmin[bb] -= one(G)
                C.block2condmax[bb] -= one(G)
            end
            for ss in (newstage + one(G)):stage
                C.stage2condmin[s] -= one(G)
            end
        end
    end
    return C
end

"""
    rowswitch_link!(newrow::G, col::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> C

Checks that `col` is currently linked and changes link from `row`, `col` to `newrow`, `col`.

Only performed if `newrow` currently unassigned and if `newrow` and `col` can be linked in the existing block.
Thus, the link block is held constant.

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function rowswitch_link!(newrow::G, col::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(C.col2row[col])
        @warn "Col not linked, no switch made"
    elseif !iszero(C.row2col[newrow])
        @warn "Newrow already assigned, no switch made"
    end
    
    stage = C.row2stage[C.col2row[col]]
    if bstages.rs2b[newrow, stage] != bstages.cs2b[col, stage]
        @warn "newrow and column not in same block in stage $stage, no switch made"
    else
        remove_link!(C.col2row[col], col, C, bstages)
        add_link!(newrow, col, stage, C, bstages)
        #C.row2col[C.col2row[col]] = zero(G)
        #C.row2block[C.col2row[col]] = zero(G)
        #C.row2stage[C.col2row[col]] = zero(G)
        #C.row2col[newrow] = col
        #C.col2row[col] = newrow
        #C.row2block[C.col2row[col]] = zero(G)
        #C.row2stage[C.col2row[col]] = zero(G)
    end
    return C
end

"""
    colswitch_link!(row::G, newcol::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> C

Checks that `row` is currently linked and changes link from `row`, `col` to `row`, `newcol`.

Only performed if `newcol` currently unassigned and if `newcol` and `row` can be linked in the existing block.
Thus, the link block is held constant.

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function colswitch_link!(row::G, newcol::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(C.row2col[row])
        @warn "Row not linked, no switch made"
    elseif !iszero(C.col2row[newcol])
        @warn "Newcol already assigned, no switch made"
    end

    stage = C.row2stage[row]
    if bstages.rs2b[row, stage] != bstages.cs2b[newcol, stage]
        @warn "row and newcolumn not in same block in stage $stage, no switch made"
    else
        C.col2row[C.row2col[row]] = zero(G)
        C.row2col[row] = newcol
        C.col2row[newcol] = row
    end
    return C
end

"""
    doubleswitch_link!(newrow::G, newcol::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> C

Performed where `newrow` currnetly linked to `col` and `row` currently linked to `newcol`.  Switches links to `newrow`, `newcol` and `row`, `col`

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function doubleswitch_link!(newrow::G, newcol::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(C.row2col[newrow])
        @warn "Newrow not linked, no switch made"
    elseif iszero(C.col2row[newcol])
        @warn "New col not linked, no switch made"
    elseif C.row2col[newrow] == newcol
        @warn "Newrow and newcol already linked, no switch made"
    end
    col = C.row2col[newrow]
    row = C.col2row[newcol]
    block = C.row2block[row]
    stage = C.row2stage[row]
    if C.row2block[row] != C.row2block[newrow]
        @warn "Existing links in different blocks, no switch made"
    elseif bstages.rs2b[newrow, stage] != block
        @warn "Newrow not in block, no switchmade"
    elseif bstages.cs2b[newcol, stage] != block
        @warn "Newcol not in block, no switchmade"
    else
        C.row2col[row] = col
        C.row2col[newrow] = newcol
        C.col2row[col] = row
        C.col2row[newcol] = newrow
    end
    return C
end

"""
    pairmove_link!(newrow::G, newcol::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> C

Determines movetype based on current link status of both `newrow` and `newcol`.

Mapping from link status to link type is as follows:
* `newrow` unlinked and `newcol` unlinked -> `add_link!`
* `newrow` unlinked and `newcol` linked -> `rowswitch_link!`
* `newrow` linked and `newcol` unlinked -> `colswitch_link!`
* `newrow` linked to `newcol` unlinked -> `remove_link!`
* `newrow` linked not to `newcol` and `newcol` linked -> `doubleswithc_link!`

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_link!`][@ref]
"""
function pairmove_link!(newrow::G, newcol::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(C.row2col[newrow])
        if iszero(C.col2row[newcol])
            return add_link!(newrow, newcol, newstage, C, bstages)
        else
            return rowswitch_link!(newrow, newcol, C, bstages)
        end
    elseif iszero(C.col2row[newcol]) #row not zero, col zero
        return  colswitch_link!(newrow, newcol, C, bstages)
    elseif C.row2col[newrow] == newcol
        return remove_link!(newrow, newcol, C, bstages)
    else
        return doubleswitch_link!(newrow, newcol, C, bstages)
    end
end

"""
    pairmove_inverse(newrow::G, newcol::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer -> revrow, revcol

Determines inverse of `pairmove_link!` so that `pairmove_link!(newrow, newcol,...)` and then `parimove_link!(revmove, revol,...)` will result in the original `C`.

See also: [`add_link!`][@ref], [`remove_link!`][@ref], [`stageswitch_link!`][@ref], [`rowswitch_link!`][@ref], [`colswitch_link!`][@ref], [`doubleswitch_link!`][@ref], [`pairmove_inverse`][@ref]
"""
function pairmove_inverse(newrow::G, newcol::G, newstage::G, C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    if iszero(C.row2col[newrow])
        if iszero(C.col2row[newcol])
            return newrow, newcol
        else
            return C.col2row[newcol], newcol
        end
    elseif iszero(C.col2row[newcol]) #row not zero, col zero
        return  newrow, C.row2col[newrow]
    elseif C.row2col[newrow] == newcol
        return newrow, newcol
    else
        return newrow, C.row2col[newrow] #could also do C.col2row[newcol], newcol
    end
end

function get_blockopenrows(block::Integer, C::StageLinkMatrix, bstages::BlockStages)
    openrowidx = iszero.(C.row2col[bstages.b2rows[block]])
    return bstages.b2rows[block][openrowidx]
end

function get_blockopencols(block::Integer, C::StageLinkMatrix, bstages::BlockStages)
    opencolidx = iszero.(C.col2row[bstages.b2cols[block]])
    return bstages.b2cols[block][opencolidx]
end

function firststage_counts(C::StageLinkMatrix{G}, bstages::BlockStages{G}) where G <: Integer
    cts = zeros(G, bstages.nstage)
    for row in 1:bstages.nrow
        if !iszero(C.row2col[row])
            cts[minimumstage_link(row, C.row2col[row], bstages)] += one(G)
        end
    end
    return cts
end

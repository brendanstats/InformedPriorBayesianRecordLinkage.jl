module InformedPriorBayesianRecordLinkage

#using StatsBase: sample
using SpecialFunctions, StatsFuns, StatsBase, SparseArrays, Random
using Distributions: Beta, Binomial, rand
using BayesianRecordLinkage: get_loglik,
    get_counts,
    lsqrt,
    lbarker,
    get_obsidxcounts,
    gibbs_MU_draw,
    ComparisonSummary,
    SparseComparisonSummary,
    ParameterChain,
    PosthocBlocks,
    sample_proposal_full,
    counts2indicies
import BayesianRecordLinkage: add_link!, #stagelinkmatrix.jl
    remove_link!, rowswitch_link!, colswitch_link!, doubleswitch_link!,
    counts_matches, #mcmc.jl
    mh_gibbs_count,
    mh_gibbs_trace

export BlockStages,
    minimumstage_link
export StageLinkMatrix,
    stageswitch_link!,
    rowswitch_link!,
    colswitch_link!,
    doubleswitch_link!,
    pairmove_link!,
    pairmove_inverse,
    get_blockopenrows,
    get_blockopencols,
    firststage_counts
export blockbetabipartite,
    lblockbetabipartite,
    blockbetabipartite_ratio,
    lblockbetabipartite_ratio,
    lblockbetabipartite_sequence,
    lblockbetabipartite_sequencep1,
    lp1,
    lm1,
    le1,
    ls1,
    lblockpriorratio,
    draw_stageblockbetabipartite
export count_blockprevlinks,
    count_stageprevlinks,
    count_prevlinks,
    pairmove_linkdelta,
    pairmove_loglikpCratio,
    pairmove_countdelta,
    pairmove_update!,
    log_move_weights,
    pairmove_locally_balanced_update!,
    stage_gibbs!,
    pairmove_globally_balanced_update!,
    pairmove_locally_balanced_sqrt_update!,
    pairmove_locally_balanced_barker_update!
export count_matches,
    restrict_blocks,
    filter_rows_stage,
    filter_cols_stage,
    mh_gibbs_count,
    mh_gibbs_trace,
    mh_gibbs_count_fixedparam

include("blockstages.jl")
include("stagelinkmatrix.jl")
include("blockbetabipartite.jl")
include("move.jl")
include("mcmc.jl")

end

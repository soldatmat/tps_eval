using BioSequences
using BioAlignments

count_positives(alignment_result::PairwiseAlignmentResult, submat::SubstitutionMatrix) = count_positives(alignment(alignment_result), submat)

function count_positives(alignment::PairwiseAlignment, submat::SubstitutionMatrix)
    seq = alignment.a.seq
    ref = alignment.b
    anchors = alignment.a.aln.anchors

    @assert anchors[1].op == OP_START
    s = anchors[1].seqpos
    a = anchors[1].alnpos
    r = anchors[1].refpos

    positive_count = 0
    for anchor = anchors[2:end]
        positive_count_addition = _count_positives(seq, ref, anchor, s, a, r, submat)
        positive_count += positive_count_addition
        s = anchor.seqpos
        a = anchor.alnpos
        r = anchor.refpos
    end

    return positive_count
end

function _count_positives(seq::T, ref::T, anchor::AlignmentAnchor, s::Int, a::Int, r::Int, submat::SubstitutionMatrix) where {T<:LongSequence}
    ismatchop(anchor.op) || return 0

    positive_count_addition = 0
    for i = 1:anchor.alnpos-a
        positive_count_addition += submat[seq[s+i], ref[r+i]] > 0
    end

    return positive_count_addition
end

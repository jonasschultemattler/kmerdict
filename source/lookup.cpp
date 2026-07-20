#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"


template<bool use_shape, bool locate>
uint64_t RSHash::lookup1(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;
    uint64_t* offsets = new uint64_t[m_thres1];
    
    uint64_t minimiser, minimiser_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & kernel_mask) >> 2*shape_overlap_right;
            kernel_rev = (kmer_rc & kernel_mask_rev) >> 2*shape_overlap_left;
            kmer = _pext_u64(kmer, shape_mask);
            kmer_rc = _pext_u64(kmer_rc, shape_mask_rev);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            if constexpr (locate) {
                if constexpr (use_shape)
                    occurences += hashmap.contains(kmer) || hashmap_rc.contains(kmer_rc);
                else
                    occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
            }
            else {
                if constexpr (use_shape)
                    occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                else
                    occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
            }
        }
    }

    delete[] offsets;

    return occurences;
}

template<bool use_shape, bool locate>
uint64_t RSHash::lookup2(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres2];
    
    uint64_t minimiser, minimiser_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & kernel_mask) >> 2*shape_overlap_right;
            kernel_rev = (kmer_rc & kernel_mask_rev) >> 2*shape_overlap_left;
            kmer = _pext_u64(kmer, shape_mask);
            kmer_rc = _pext_u64(kmer_rc, shape_mask_rev);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                occurences += check<2, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else {
                if constexpr (locate) {
                    if constexpr (use_shape)
                        occurences += hashmap.contains(kmer) || hashmap_rc.contains(kmer_rc);
                    else
                        occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
                }
                else {
                    if constexpr (use_shape)
                        occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                    else
                        occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
                }
            }
        }
    }

    delete[] offsets;

    return occurences;
}

template<bool use_shape, bool locate>
uint64_t RSHash::lookup3(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres3];
    
    uint64_t minimiser, minimiser_rank;
    uint64_t kernel, kernel_rev;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        uint64_t kmer_rc = crc(kmer, window_size);

        if constexpr (use_shape) {
            kernel = (kmer & kernel_mask) >> 2*shape_overlap_right;
            kernel_rev = (kmer_rc & kernel_mask_rev) >> 2*shape_overlap_left;
            kmer = _pext_u64(kmer, shape_mask);
            kmer_rc = _pext_u64(kmer_rc, shape_mask_rev);
        }
        else {
            kernel = kmer;
            kernel_rev = kmer_rc;
        }

        minimiser = find_minimiser<1>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);

        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                occurences += check<2, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else {
                minimiser = find_minimiser<3>(kernel, kernel_rev, left_minimiser_position, right_minimiser_position);
                if(r3.contains(minimiser, minimiser_rank)) {
                    size_t p = s3_select.select(minimiser_rank);
                    size_t no_minimiser = s3_select.select(minimiser_rank+1) - p;
                    occurences += check<3, use_shape>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
                }
                else {
                    if constexpr (locate) {
                        if constexpr (use_shape)
                            occurences += hashmap.contains(kmer) || hashmap_rc.contains(kmer_rc);
                        else
                            occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
                    }
                    else {
                        if constexpr (use_shape)
                            occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                        else
                            occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
                    }
                }
            }
        }
    }

    delete[] offsets;

    return occurences;
}


uint64_t RSHash::lookup(const std::vector<uint64_t> &kmers)
{
    bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();
    if(level == 1)
        if(use_shape)
            if(loc)
                return lookup1<true, true>(kmers);
            else
                return lookup1<true, false>(kmers);
        else
            if(loc)
                return lookup1<false, true>(kmers);
            else
                return lookup1<false, false>(kmers);
    else if(level == 2)
        if(use_shape)
            if(loc)
                return lookup2<true, true>(kmers);
            else
                return lookup2<true, false>(kmers);
        else
            if(loc)
                return lookup2<false, true>(kmers);
            else
                return lookup2<false, false>(kmers);
    else if(level == 3)
        if(use_shape)
            if(loc)
                return lookup3<true, true>(kmers);
            else
                return lookup3<true, false>(kmers);
        else
            if(loc)
                return lookup3<false, true>(kmers);
            else
                return lookup3<false, false>(kmers);
    else
        return 0;
}


template<int level, bool use_shape>
inline bool RSHash::check(const uint64_t kmer, const uint64_t kmer_rc,
    uint64_t* offsets, const size_t p, const size_t no_skmers,
    const size_t left_minimiser_position, const size_t right_minimiser_position)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;
    
    for(size_t i = 0; i < no_skmers; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i)-span+1;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i)-span+1;
        if constexpr (level == 3)
            offsets[i] = offsets3.access(p+i)-span+1;
    }

    for(size_t i = 0; i < no_skmers; i++) {
        uint64_t offset = offsets[i];
        uint64_t pos = span-1-left_minimiser_position;
        uint64_t pos_rc = left_minimiser_position;

        uint64_t hash_fwd = get_word64(offset + pos - shape_overlap_right) & windowmask;
        uint64_t hash_rc = get_word64(offset + pos_rc - shape_overlap_left) & windowmask;

        if constexpr (use_shape) {
            hash_fwd = _pext_u64(hash_fwd, shape_mask);
            hash_rc = _pext_u64(hash_rc, shape_mask_rev);
        }

        if(kmer == hash_fwd || kmer_rc == hash_rc)
            return true;

        if(left_minimiser_position != k-m1-right_minimiser_position) {
            pos = right_minimiser_position;
            pos_rc = span-1-right_minimiser_position;

            hash_fwd = get_word64(offset + pos - shape_overlap_right) & windowmask;
            hash_rc = get_word64(offset + pos_rc - shape_overlap_left) & windowmask;

            if constexpr (use_shape) {
                hash_fwd = _pext_u64(hash_fwd, shape_mask);
                hash_rc = _pext_u64(hash_rc, shape_mask_rev);
            }

            if(kmer == hash_fwd || kmer_rc == hash_rc)
                return true;
        }

    }

    return false;
}

template<bool use_shape>
inline bool RSHash::extend_in_text(uint64_t &text_pos, uint64_t start, uint64_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, uint64_t &window, uint64_t &window_rev)
{
    if(forward) {
        if(++text_pos < end) {
            const uint64_t new_rank = get_base(text_pos);
            if constexpr (use_shape) {
                window = (window >> 2) | (new_rank << windowshift);
                return _pext_u64(window, shape_mask) == query;
            }
            else
                return new_rank == (query >> windowshift);
        }
    }
    else {
        if(--text_pos >= start) {
            const uint64_t new_rank = get_base(text_pos);
            if constexpr (use_shape) {
                window_rev = ((window_rev << 2) | new_rank);
                return _pext_u64(window_rev, shape_mask_rev) == query_rc;
            }
            else
                return new_rank == (query_rc & 0b11);
        }
    }
    return false;
}


template<int level, bool use_shape>
inline bool RSHash::check_minimiser_pos(uint64_t *buffer, uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos,
    uint64_t &text_kmer, uint64_t &text_kmer_rc)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;
    
    uint64_t candidate, candidate_rc, window, window_rev, pos, pos_rc;
    if constexpr (use_shape) {
        pos = overlap - shape_overlap_right + span-1-minimiser_pos;
        pos_rc = overlap - shape_overlap_left + minimiser_pos;
        offset -= overlap;
        window = buffer[s + pos];
        window_rev = buffer[s + pos_rc];
        candidate = _pext_u64(window, shape_mask);
        candidate_rc = _pext_u64(window_rev, shape_mask_rev);
    }
    else {
        pos = span-1-minimiser_pos;
        pos_rc = minimiser_pos;
        window = candidate = buffer[s+pos];
        window_rev = candidate_rc = buffer[s+pos_rc];
    }

    if(candidate == query && check_overlap(offset + pos, start_pos, end_pos)) {
        forward = true;
        text_pos = offset + pos + window_size - 1;
        text_kmer = window;
        return true;
    }
    if(candidate_rc == queryrc && check_overlap(offset + pos_rc, start_pos, end_pos)) {
        forward = false;
        text_pos = offset + pos_rc;
        text_kmer_rc = window_rev;
        return true;
    }

    return false;
}


template<int level, bool use_shape>
inline bool RSHash::check_minimiser_pos2(uint64_t *buffer, uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos,
    uint64_t &text_kmer, uint64_t &text_kmer_rc)
{   
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;
    
    uint64_t left_candidate, left_candidate_rc, right_candidate, right_candidate_rc;
    uint64_t left_window, left_window_rev, right_window, right_window_rev;
    uint64_t left_pos, left_pos_rc, right_pos, right_pos_rc;
    if constexpr (use_shape) {
        left_pos = overlap - shape_overlap_right + span-1-left_minimiser_pos;
        left_pos_rc = overlap - shape_overlap_left + left_minimiser_pos;
        right_pos = overlap - shape_overlap_right + right_minimiser_pos;
        right_pos_rc = overlap - shape_overlap_left + span-1-right_minimiser_pos;
        left_window = buffer[s + left_pos];
        left_window_rev = buffer[s + left_pos_rc];
        right_window = buffer[s + right_pos];
        right_window_rev = buffer[s + right_pos_rc];
        left_candidate = _pext_u64(left_window, shape_mask);
        left_candidate_rc = _pext_u64(left_window_rev, shape_mask_rev);
        right_candidate = _pext_u64(right_window, shape_mask);
        right_candidate_rc = _pext_u64(right_window_rev, shape_mask_rev);
        offset -= overlap;
    }
    else {
        left_pos = span-1-left_minimiser_pos;
        left_pos_rc = left_minimiser_pos;
        right_pos = right_minimiser_pos;
        right_pos_rc = span-1-right_minimiser_pos;
        left_window = left_candidate = buffer[s + left_pos];
        left_window_rev = left_candidate_rc = buffer[s + left_pos_rc];
        right_window = right_candidate = buffer[s + right_pos];
        right_window_rev = right_candidate_rc = buffer[s + right_pos_rc];
    }

    if(left_candidate == query && check_overlap(offset + left_pos, start_pos, end_pos)) {
        forward = true;
        text_pos = offset + left_pos + window_size - 1;
        text_kmer = left_window;
        return true;
    }
    if(left_candidate_rc == queryrc && check_overlap(offset + left_pos_rc, start_pos, end_pos)) {
        forward = false;
        text_pos = offset + left_pos_rc;
        text_kmer_rc = left_window_rev;
        return true;
    }
    if(right_candidate == query && check_overlap(offset + right_pos, start_pos, end_pos)) {
        forward = true;
        text_pos = offset + right_pos + window_size - 1;
        text_kmer = right_window;
        return true;
    }
    if(right_candidate_rc == queryrc && check_overlap(offset + right_pos_rc, start_pos, end_pos)) {
        forward = false;
        text_pos = offset + right_pos_rc;
        text_kmer = right_window_rev;
        return true;
    }

    return false;
}


template<int level, bool use_shape>
inline bool RSHash::lookup_buffer(uint64_t* buffer, uint64_t *offsets, const size_t no_skmers,
    const uint64_t query, const uint64_t queryrc,
    uint64_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &start_pos, uint64_t &end_pos, uint64_t &text_kmer, uint64_t &text_kmer_rc)
{
    size_t span, m;
    if constexpr (level == 1) {
        span = span1;
        m = m1;
    }
    if constexpr (level == 2) {
        span = span2;
        m = m2;
    }
    if constexpr (level == 3) {
        span = span3;
        m = m3;
    }

    size_t s = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_skmers; i++) {
            if(check_minimiser_pos2<level, use_shape>(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos, text_kmer, text_kmer_rc))
                return true;
            s += span;
            if constexpr (use_shape)
                s += overlap;
        }
    }
    else {
        for(size_t i = 0; i < no_skmers; i++) {
            if(check_minimiser_pos<level, use_shape>(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, forward, text_pos, start_pos, end_pos, text_kmer, text_kmer_rc))
                return true;
            s += span;
            if constexpr (use_shape)
                s += overlap;
        }
    }
    
    return false;
}


uint64_t RSHash::streaming_lookup(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();
    if(level == 1) {
        if(use_shape)
            if(loc)
                return streaming_lookup1<true, true>(query, extensions);
            else
                return streaming_lookup1<true, false>(query, extensions);
        else
            if(loc)
                return streaming_lookup1<false, true>(query, extensions);
            else
                return streaming_lookup1<false, false>(query, extensions);
    }
    else if(level == 2) {
        if(use_shape)
            if(loc)
                return streaming_lookup2<true, true>(query, extensions);
            else
                return streaming_lookup2<true, false>(query, extensions);
        else
            if(loc)
                return streaming_lookup2<false, true>(query, extensions);
            else
                return streaming_lookup2<false, false>(query, extensions);
    }
    else if(level == 3) {
        if(use_shape)
            if(loc)
                return streaming_lookup3<true, true>(query, extensions);
            else
                return streaming_lookup3<true, false>(query, extensions);
        else
            if(loc)
                return streaming_lookup3<false, true>(query, extensions);
            else
                return streaming_lookup3<false, false>(query, extensions);
    }
    else
        return 0;
}


template<bool use_shape, bool locate>
uint64_t RSHash::streaming_lookup1(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF;
    uint64_t current_neg_minimiser1=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* buffer1 = new uint64_t[m_thres1 * (span1 + overlap)];
    size_t no_skmers1;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    uint64_t kmer, kmer_rc, kernel, kernel_rev;
    uint64_t text_kmer, text_kmer_rc;

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kmer = _pext_u64(window.value, shape_mask);
            kmer_rc = _pext_u64(window.value_rev, shape_mask_rev);
        }
        else {
            kmer = window.value;
            kmer_rc = window.value_rev;
        }

        if(found && extend_in_text<use_shape>(text_pos, sequence_begin, sequence_end, forward, kmer, kmer_rc, text_kmer, text_kmer_rc)) {
            occurences++;
            extensions++;
            rolling = false;
        }
        else {
            if constexpr (use_shape) {
                kernel = (window.value & kernel_mask) >> 2*shape_overlap_right;
                kernel_rev = (window.value_rev & kernel_mask_rev) >> 2*shape_overlap_left;
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, use_shape>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
                current_minimiser1 = minimiser1;
            }    
            else {
                if constexpr (use_shape)
                    occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                else
                    occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
                
                found = false;
                current_neg_minimiser1 = minimiser1;
            }
        }
    }

    delete[] offsets1;
    delete[] buffer1;
    
    return occurences;
}

template<bool use_shape, bool locate>
uint64_t RSHash::streaming_lookup2(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* buffer1 = new uint64_t[m_thres1 * (span1 + overlap)];
    uint64_t* buffer2 = new uint64_t[m_thres2 * (span2 + overlap)];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;
    uint64_t kmer, kmer_rc, kernel, kernel_rev;
    uint64_t text_kmer, text_kmer_rc;

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kmer = _pext_u64(window.value, shape_mask);
            kmer_rc = _pext_u64(window.value_rev, shape_mask_rev);
        }
        else {
            kmer = window.value;
            kmer_rc = window.value_rev;
        }

        if(found && extend_in_text<use_shape>(text_pos, sequence_begin, sequence_end, forward, kmer, kmer_rc, text_kmer, text_kmer_rc)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
        }
        else {
            if constexpr (use_shape) {
                kernel = (window.value & kernel_mask) >> 2*shape_overlap_right;
                kernel_rev = (window.value_rev & kernel_mask_rev) >> 2*shape_overlap_left;
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling1)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, use_shape>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer<2, use_shape>(buffer2, offsets2, no_skmers2, kmer, kmer_rc, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                    occurences += found;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer<2, use_shape>(offsets2, buffer2, p, no_skmers2);
                    found = lookup_buffer<2, use_shape>(buffer2, offsets2, no_skmers2, kmer, kmer_rc, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }   
                else {
                    if constexpr (use_shape) // todo: and shape is not symmetric
                        occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                    else
                        occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
                    found = false;
                    current_neg_minimiser1 = minimiser1;
                    current_neg_minimiser2 = minimiser2;
                }
            }
        }
    }

    delete[] offsets1;
    delete[] offsets2;
    delete[] buffer1;
    delete[] buffer2;
    
    return occurences;
}

template<bool use_shape, bool locate>
uint64_t RSHash::streaming_lookup3(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF, current_minimiser3=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF, current_neg_minimiser3=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* buffer1 = new uint64_t[m_thres1 * (span1 + overlap)];
    uint64_t* buffer2 = new uint64_t[m_thres2 * (span2 + overlap)];
    uint64_t* buffer3 = new uint64_t[m_thres3 * (span3 + overlap)];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t sequence_begin, sequence_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    bool rolling3 = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;
    size_t left_minimiser3_position, right_minimiser3_position;
    uint64_t minimiser3, minimiser3_rank;
    uint64_t text_kmer, text_kmer_rc;
    uint64_t kmer, kmer_rc, kernel, kernel_rev;

    uint64_t occurences = 0;
    for(auto && window : query | rshash::views::kmerview({.window_size = window_size}))
    {
        if constexpr (use_shape) {
            kmer = _pext_u64(window.value, shape_mask);
            kmer_rc = _pext_u64(window.value_rev, shape_mask_rev);
        }
        else {
            kmer = window.value;
            kmer_rc = window.value_rev;
        }

        if(found && extend_in_text<use_shape>(text_pos, sequence_begin, sequence_end, forward, kmer, kmer_rc, text_kmer, text_kmer_rc)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
            rolling3 = false;
        }
        else {
            if constexpr (use_shape) {
                kernel = (window.value & kernel_mask) >> 2*shape_overlap_right;
                kernel_rev = (window.value_rev & kernel_mask_rev) >> 2*shape_overlap_left;
            }
            else {
                kernel = window.value;
                kernel_rev = window.value_rev;
            }

            if(rolling1)
                update_minimiser<1>(kernel, kernel_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kernel, kernel_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
                rolling2 = false;
                rolling3 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, use_shape>(offsets1, buffer1, p, no_skmers1);
                found = lookup_buffer<1, use_shape>(buffer1, offsets1, no_skmers1, kmer, kmer_rc, text_pos, left_minimiser1_position, right_minimiser1_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
                rolling3 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kernel, kernel_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kernel, kernel_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer<2, use_shape>(buffer2, offsets2, no_skmers2, kmer, kmer_rc, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                    occurences += found;
                    rolling3 = false;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer<2, use_shape>(offsets2, buffer2, p, no_skmers2);
                    found = lookup_buffer<2, use_shape>(buffer2, offsets2, no_skmers2, kmer, kmer_rc, text_pos, left_minimiser2_position, right_minimiser2_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                    rolling3 = false;
                }
                else {
                    if(rolling3)
                        update_minimiser<3>(kernel, kernel_rev, minimiser3, left_minimiser3_position, right_minimiser3_position);
                    else {
                        minimiser3 = find_minimiser<3>(kernel, kernel_rev, left_minimiser3_position, right_minimiser3_position);
                        rolling3 = true;
                    }

                    if(minimiser3 == current_minimiser3) {
                        found = lookup_buffer<3, use_shape>(buffer3, offsets3, no_skmers3, kmer, kmer_rc, text_pos, left_minimiser3_position, right_minimiser3_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                        occurences += found;
                    }
                    else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser3_rank)) {
                        const size_t p = s3_select.select(minimiser3_rank);
                        no_skmers3 = s3_select.select(minimiser3_rank+1) - p;

                        fill_buffer<3, use_shape>(offsets3, buffer3, p, no_skmers3);
                        found = lookup_buffer<3, use_shape>(buffer3, offsets3, no_skmers3, kmer, kmer_rc, text_pos, left_minimiser3_position, right_minimiser3_position, forward, sequence_begin, sequence_end, text_kmer, text_kmer_rc);
                        occurences += found;
                        current_minimiser3 = minimiser3;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                    }
                    else {
                        if constexpr (use_shape)
                            occurences += hashset.contains(kmer) || hashset_rc.contains(kmer_rc);
                        else
                            occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
                        found = false;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                        current_neg_minimiser3 = minimiser3;
                    }
                }
            }
        }
    }

    delete[] offsets1;
    delete[] offsets2;
    delete[] offsets3;
    delete[] buffer1;
    delete[] buffer2;
    delete[] buffer3;
    
    return occurences;
}

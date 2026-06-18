#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"


uint64_t RSHash::lookup1(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;
    uint64_t* offsets = new uint64_t[m_thres1];
    
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        const uint64_t kmer_rc = crc(kmer, k);
        minimiser = find_minimiser<1>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else
            occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
    }

    delete[] offsets;

    return occurences;
}

uint64_t RSHash::lookup2(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres2];
    
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        const uint64_t kmer_rc = crc(kmer, k);
        minimiser = find_minimiser<1>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                occurences += check<2>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else
                occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
        }
    }

    delete[] offsets;

    return occurences;
}

uint64_t RSHash::lookup3(const std::vector<uint64_t> &kmers)
{
    uint64_t occurences = 0;

    uint64_t* offsets = new uint64_t[m_thres3];
    
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(uint64_t kmer : kmers) {
        const uint64_t kmer_rc = crc(kmer, k);
        minimiser = find_minimiser<1>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);

        if(r1.contains(minimiser, minimiser_rank)) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_minimiser = s1_select.select(minimiser_rank+1) - p;
            occurences += check<1>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
        }
        else {
            minimiser = find_minimiser<2>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
            if(r2.contains(minimiser, minimiser_rank)) {
                size_t p = s2_select.select(minimiser_rank);
                size_t no_minimiser = s2_select.select(minimiser_rank+1) - p;
                occurences += check<2>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
            }
            else {
                minimiser = find_minimiser<3>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
                if(r3.contains(minimiser, minimiser_rank)) {
                    size_t p = s3_select.select(minimiser_rank);
                    size_t no_minimiser = s3_select.select(minimiser_rank+1) - p;
                    occurences += check<3>(kmer, kmer_rc, offsets, p, no_minimiser, left_minimiser_position, right_minimiser_position);
                }
                else
                    occurences += hashset.contains(std::min<uint64_t>(kmer, kmer_rc));
            }
        }
    }

    delete[] offsets;

    return occurences;
}

uint64_t RSHash::lookup(const std::vector<uint64_t> &kmers)
{
    if(level == 1)
        return lookup1(kmers);
    else if(level == 2)
        return lookup2(kmers);
    else if(level == 3)
        return lookup3(kmers);
    else
        return 0;
}


template<int level>
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
        const uint64_t o = offsets[i];

        uint64_t hash_rc = get_word64(o + left_minimiser_position) & kmermask;
        uint64_t hash_fwd = get_word64(o + span-1-left_minimiser_position) & kmermask;

        if(kmer == hash_fwd || kmer_rc == hash_rc)
            return true;

        if(left_minimiser_position != k-m1-right_minimiser_position) {
            hash_fwd = get_word64(o + right_minimiser_position) & kmermask;
            hash_rc = get_word64(o + span-1-right_minimiser_position) & kmermask;

            if(kmer == hash_fwd || kmer_rc == hash_rc)
                return true;
        }

    }

    return false;
}


inline bool RSHash::extend_in_text(uint64_t &text_pos, uint64_t start, uint64_t end,
    bool forward, const uint64_t query, const uint64_t query_rc)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query >> (2*(k-1))));
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query_rc & 0b11));
            return found;
        }
    }
    return false;
}


template<int level>
inline bool RSHash::check_minimiser_pos(uint64_t *buffer, const uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t e, const size_t minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-minimiser_pos] == query) {
        forward = true;
        text_pos = offset + e-1-s-minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }

    return false;
}


template<int level>
inline bool RSHash::check_minimiser_pos2(uint64_t *buffer, const uint64_t offset,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t e, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    if(buffer[s+left_minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + left_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-left_minimiser_pos] == query) {
        forward = true;
        text_pos = offset + e-1-s-left_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+right_minimiser_pos] == query) {
        forward = true;
        text_pos = offset + right_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[e-1-right_minimiser_pos] == queryrc) {
        forward = false;
        text_pos = offset + e-1-s-right_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos))
            return true;
    }

    return false;
}


template<int level>
inline bool RSHash::lookup_buffer(uint64_t* buffer, uint64_t *offsets, const size_t no_skmers,
    const uint64_t query, const uint64_t queryrc,
    uint64_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, uint64_t &start_pos, uint64_t &end_pos)
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

    size_t s = 0, e = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;
            if(check_minimiser_pos2<level>(buffer, offsets[i], query, queryrc, s, e, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    else {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;
            if(check_minimiser_pos<level>(buffer, offsets[i], query, queryrc, s, e, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    
    return false;
}

uint64_t RSHash::streaming_lookup(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    if(level == 1)
        return streaming_lookup1(query, extensions);
    else if(level == 2)
        return streaming_lookup2(query, extensions);
    else if(level == 3)
        return streaming_lookup3(query, extensions);
    else
        return 0;
}


uint64_t RSHash::streaming_lookup1(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    const uint64_t shift = 2*(k-1);
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF;
    uint64_t current_neg_minimiser1=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* buffer1 = new uint64_t[m_thres1 * span1];
    size_t no_skmers1;
    uint64_t unitig_begin, unitig_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;

    uint64_t occurences = 0;
    for(auto && kmer : query | rshash::views::kmerview({.window_size = k}))
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, kmer.value, kmer.value_rev)) {
            occurences++;
            extensions++;
            rolling = false;
        }
        else {
            if(rolling)
                update_minimiser<1>(kmer.value, kmer.value_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kmer.value, kmer.value_rev, left_minimiser1_position, right_minimiser1_position);
                rolling = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, false>(offsets1, buffer1, p, no_skmers1, shift);
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                current_minimiser1 = minimiser1;
            }    
            else {
                occurences += hashset.contains(std::min<uint64_t>(kmer.value, kmer.value_rev));
                found = false;
                current_neg_minimiser1 = minimiser1;
            }
        }
    }

    delete[] offsets1;
    delete[] buffer1;
    
    return occurences;
}


uint64_t RSHash::streaming_lookup2(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = rshash::views::kmerview({.window_size = k});

    const uint64_t shift = 2*(k-1);
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* buffer1 = new uint64_t[m_thres1 * span1];
    uint64_t* buffer2 = new uint64_t[m_thres2 * span2];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t unitig_begin, unitig_end;
    uint64_t text_pos;
    bool forward;
    bool found = false;
    bool rolling1 = false;
    bool rolling2 = false;
    size_t left_minimiser1_position, right_minimiser1_position;
    uint64_t minimiser1, minimiser1_rank;
    size_t left_minimiser2_position, right_minimiser2_position;
    uint64_t minimiser2, minimiser2_rank;
    uint64_t occurences = 0;

    for(auto && kmer : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, kmer.value, kmer.value_rev)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
        }
        else {
            if(rolling1)
                update_minimiser<1>(kmer.value, kmer.value_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kmer.value, kmer.value_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, false>(offsets1, buffer1, p, no_skmers1, shift);
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kmer.value, kmer.value_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kmer.value, kmer.value_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, kmer.value, kmer.value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer<2, false>(offsets2, buffer2, p, no_skmers2, shift);
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, kmer.value, kmer.value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }   
                else {
                    occurences += hashset.contains(std::min<uint64_t>(kmer.value, kmer.value_rev));
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


uint64_t RSHash::streaming_lookup3(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = rshash::views::kmerview({.window_size = k});

    const uint64_t shift = 2*(k-1);
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF, current_minimiser3=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF, current_neg_minimiser3=INF;
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* buffer1 = new uint64_t[m_thres1 * span1];
    uint64_t* buffer2 = new uint64_t[m_thres2 * span2];
    uint64_t* buffer3 = new uint64_t[m_thres3 * span3];
    size_t no_skmers1, no_skmers2, no_skmers3;
    uint64_t unitig_begin, unitig_end;
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
    uint64_t occurences = 0;

    for(auto && kmer : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, kmer.value, kmer.value_rev)) {
            occurences++;
            extensions++;
            rolling1 = false;
            rolling2 = false;
            rolling3 = false;
        }
        else {
            if(rolling1)
                update_minimiser<1>(kmer.value, kmer.value_rev, minimiser1, left_minimiser1_position, right_minimiser1_position);
            else {
                minimiser1 = find_minimiser<1>(kmer.value, kmer.value_rev, left_minimiser1_position, right_minimiser1_position);
                rolling1 = true;
            }

            if(minimiser1 == current_minimiser1) {
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                rolling2 = false;
                rolling3 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser1_rank)) {
                const size_t p = s1_select.select(minimiser1_rank);
                no_skmers1 = s1_select.select(minimiser1_rank+1) - p;

                fill_buffer<1, false>(offsets1, buffer1, p, no_skmers1, shift);
                found = lookup_buffer<1>(buffer1, offsets1, no_skmers1, kmer.value, kmer.value_rev, text_pos, left_minimiser1_position, right_minimiser1_position, forward, unitig_begin, unitig_end);
                occurences += found;
                current_minimiser1 = minimiser1;
                rolling2 = false;
                rolling3 = false;
            }
            else {
                if(rolling2)
                    update_minimiser<2>(kmer.value, kmer.value_rev, minimiser2, left_minimiser2_position, right_minimiser2_position);
                else {
                    minimiser2 = find_minimiser<2>(kmer.value, kmer.value_rev, left_minimiser2_position, right_minimiser2_position);
                    rolling2 = true;
                }

                if(minimiser2 == current_minimiser2) {
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, kmer.value, kmer.value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                    rolling3 = false;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser2_rank)) {
                    const size_t p = s2_select.select(minimiser2_rank);
                    no_skmers2 = s2_select.select(minimiser2_rank+1) - p;

                    fill_buffer<2, false>(offsets2, buffer2, p, no_skmers2, shift);
                    found = lookup_buffer<2>(buffer2, offsets2, no_skmers2, kmer.value, kmer.value_rev, text_pos, left_minimiser2_position, right_minimiser2_position, forward, unitig_begin, unitig_end);
                    occurences += found;
                    current_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                    rolling3 = false;
                }
                else {
                    if(rolling3)
                        update_minimiser<3>(kmer.value, kmer.value_rev, minimiser3, left_minimiser3_position, right_minimiser3_position);
                    else {
                        minimiser3 = find_minimiser<3>(kmer.value, kmer.value_rev, left_minimiser3_position, right_minimiser3_position);
                        rolling3 = true;
                    }

                    if(minimiser3 == current_minimiser3) {
                        found = lookup_buffer<3>(buffer3, offsets3, no_skmers3, kmer.value, kmer.value_rev, text_pos, left_minimiser3_position, right_minimiser3_position, forward, unitig_begin, unitig_end);
                        occurences += found;
                    }
                    else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser3_rank)) {
                        const size_t p = s3_select.select(minimiser3_rank);
                        no_skmers3 = s3_select.select(minimiser3_rank+1) - p;

                        fill_buffer<3, false>(offsets3, buffer3, p, no_skmers3, shift);
                        found = lookup_buffer<3>(buffer3, offsets3, no_skmers3, kmer.value, kmer.value_rev, text_pos, left_minimiser3_position, right_minimiser3_position, forward, unitig_begin, unitig_end);
                        occurences += found;
                        current_minimiser3 = minimiser3;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                    }
                    else {
                        occurences += hashset.contains(std::min<uint64_t>(kmer.value, kmer.value_rev));
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


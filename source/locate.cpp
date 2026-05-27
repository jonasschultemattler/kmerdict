#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>
#include "rshash.hpp"


template<int level>
inline bool RSHash::report_minimiser_pos(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, size_t &i, const size_t s,
    const size_t minimiser_pos,
    std::vector<std::pair<uint64_t, bool>> &positions)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    bool found = false;
    uint64_t start_pos, end_pos;
    if(buffer[s+minimiser_pos] == kmerrc) {
        const uint64_t text_pos = offset + minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos)) {
            // positions[i++] = {text_pos, true};
            i++;
            found = true;
        }
    }
    if(buffer[s+span-1-minimiser_pos] == kmer) {
        const uint64_t text_pos = offset + span-1-minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos)) {
            // positions[i++] = {text_pos, false};
            i++;
            found = true;
        }
    }
    return found;
}

template<int level>
inline bool RSHash::report_minimiser_pos2(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, size_t &i, const size_t s,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<std::pair<uint64_t, bool>> &positions)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    else if constexpr (level == 2)
        span = span2;
    else if constexpr (level == 3)
        span = span3;

    bool found = false;
    uint64_t start_pos, end_pos;
    
    if(buffer[s+left_minimiser_pos] == kmerrc) {
        const uint64_t text_pos = offset + left_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos)) {
            // positions[i++] = {text_pos, true};
            i++;
            found = true;
        }
    }
    if(buffer[s+span-1-left_minimiser_pos] == kmer) {
        const uint64_t text_pos = offset + span-1-left_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos)) {
            // positions[i++] = {text_pos, false};
            i++;
            found = true;
        }
    }
    if(buffer[s+right_minimiser_pos] == kmer) {
        const uint64_t text_pos = offset + right_minimiser_pos + k - 1;
        if(check_overlap<level>(offset, text_pos-k+1, start_pos, end_pos)) {
            // positions[i++] = {text_pos, true};
            i++;
            found = true;
        }
    }
    if(buffer[s+span-1-right_minimiser_pos] == kmerrc) {
        const uint64_t text_pos = offset + span-1-right_minimiser_pos;
        if(check_overlap<level>(offset, text_pos, start_pos, end_pos)) {
            // positions[i++] = {text_pos, false};
            i++;
            found = true;
        }
    }
    return found;
}

template<int level>
inline void RSHash::locate_buffer(uint64_t *buffer, uint64_t *offsets, const size_t no_minimiser,
    const uint64_t query, const uint64_t queryrc,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    std::vector<std::pair<uint64_t, bool>> &positions, size_t &found_positions, size_t &found_kmers)
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

    // todo: SIMD buffer check
    bool found = false;
    size_t s = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_minimiser; i++) {
            found |= report_minimiser_pos2<level>(buffer, offsets[i], query, queryrc, found_positions, s, left_minimiser_pos, right_minimiser_pos, positions);
            s += span;
        }
    }
    else {
        for(size_t i = 0; i < no_minimiser; i++) {
            found |= report_minimiser_pos<level>(buffer, offsets[i], query, queryrc, found_positions, s, left_minimiser_pos, positions);
            s += span;
        }
    }
    found_kmers += found;
}

size_t RSHash::streaming_locate(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<std::pair<uint64_t, bool>> &positions, size_t &found_positions)
{
    if(level == 1)
        return streaming_locate1(query, positions, found_positions);
    else if(level == 2)
        return streaming_locate2(query, positions, found_positions);
    else if(level == 3)
        return streaming_locate3(query, positions, found_positions);
    else
        return 0;
}

size_t RSHash::streaming_locate1(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<std::pair<uint64_t, bool>> &positions, size_t &found_positions)
{
    size_t found_kmers = 0;
    uint64_t current_pos_minimiser=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser=std::numeric_limits<uint64_t>::max();
    const uint64_t shift = 2*(k-1);
    uint64_t* offsets = new uint64_t[m_thres1];
    uint64_t* kmer_buffer = new uint64_t[(m_thres1) * span1];
    size_t no_minimiser;
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;
    bool begin = true;

    auto process = [&](auto&& rng)
    {
        for (auto&& kmer : rng) {
            if(begin) {
                minimiser = find_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position, right_minimiser_position);
                begin = false;
            }
            else
                update_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, minimiser, left_minimiser_position, right_minimiser_position);

            if(minimiser == current_pos_minimiser) {
                locate_buffer<1>(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position, right_minimiser_position, positions, found_positions, found_kmers);
            }
            else if(minimiser != current_neg_minimiser && r1.contains(minimiser, minimiser_rank)) {
                const size_t minimiser_position = s1_select.select(minimiser_rank);
                no_minimiser = s1_select.select(minimiser_rank+1) - minimiser_position;

                fill_buffer<1>(offsets, kmer_buffer, minimiser_position, no_minimiser, shift);
                locate_buffer<1>(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position, right_minimiser_position, positions, found_positions, found_kmers);
                current_pos_minimiser = minimiser;
            }
            else {
                if (auto it = hashmap.find(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev)); it != hashmap.end()) {
                    for(auto pos : it->second)
                        found_positions++;
                    found_kmers++;
                }
                current_neg_minimiser = minimiser;
            }
        }
    };

    if (shape != std::numeric_limits<uint32_t>::max())
        process(query | rshash::views::shapeview({.shape = shape}));
    else
        process(query | rshash::views::kmerview({.window_size = k}));


    delete[] kmer_buffer;
    delete[] offsets;

    return found_kmers;
}


size_t RSHash::streaming_locate2(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<std::pair<uint64_t, bool>> &positions, size_t &found_positions)
{
    size_t found_kmers = 0;
    uint64_t current_pos_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser2=std::numeric_limits<uint64_t>::max();
    const uint64_t shift = 2*(k-1);
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* kmer_buffer1 = new uint64_t[(m_thres1) * span1];
    uint64_t* kmer_buffer2 = new uint64_t[(m_thres2) * span2];
    size_t no_minimiser1, no_minimiser2;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    bool rolling1 = false;
    bool rolling2 = false;

    auto process = [&](auto&& rng)
    {
        for (auto&& kmer : rng) {
            if(!rolling1) {
                minimiser1 = find_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1);
                rolling1 = true;
            }
            else
                update_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

            if(minimiser1 == current_pos_minimiser1) {
                locate_buffer<1>(kmer_buffer1, offsets1, no_minimiser1, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1, positions, found_positions, found_kmers);
                rolling2 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
                const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
                no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

                fill_buffer<1>(offsets1, kmer_buffer1, minimiser_position1, no_minimiser1, shift);
                locate_buffer<1>(kmer_buffer1, offsets1, no_minimiser1, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1, positions, found_positions, found_kmers);
                current_pos_minimiser1 = minimiser1;
                rolling2 = false;
            }
            else {
                if(!rolling2) {
                    minimiser2 = find_minimiser<2>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2);
                    rolling2 = true;
                }
                else
                    update_minimiser<2>(kmer.kmer_value, kmer.kmer_value_rev, minimiser2, left_minimiser_position2, right_minimiser_position2);

                if (minimiser2 == current_pos_minimiser2) {
                    locate_buffer<2>(kmer_buffer2, offsets2, no_minimiser2, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2, positions, found_positions, found_kmers);
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                    const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                    no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                    fill_buffer<2>(offsets2, kmer_buffer2, minimiser_position2, no_minimiser2, shift);
                    locate_buffer<2>(kmer_buffer2, offsets2, no_minimiser2, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2, positions, found_positions, found_kmers);
                    current_pos_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                }
                else {
                    if (auto it = hashmap.find(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev)); it != hashmap.end()) {
                        for(auto pos : it->second)
                            found_positions++;
                        found_kmers++;
                    }
                    current_neg_minimiser1 = minimiser1;
                    current_neg_minimiser2 = minimiser2;
                }
            }
        }
    };

    if(shape != std::numeric_limits<uint32_t>::max())
        process(query | rshash::views::shapeview({.shape = shape}));
    else
        process(query | rshash::views::kmerview({.window_size = k}));

    delete[] kmer_buffer1;
    delete[] kmer_buffer2;
    delete[] offsets1;
    delete[] offsets2;

    return found_kmers;
}


size_t RSHash::streaming_locate3(const seqan3::bitpacked_sequence<seqan3::dna4> &query,
    std::vector<std::pair<uint64_t, bool>> &positions, size_t &found_positions)
{
    size_t found_kmers = 0;
    uint64_t current_pos_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_pos_minimiser3=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser3=std::numeric_limits<uint64_t>::max();
    const uint64_t shift = 2*(k-1);
    uint64_t* offsets1 = new uint64_t[m_thres1];
    uint64_t* offsets2 = new uint64_t[m_thres2];
    uint64_t* offsets3 = new uint64_t[m_thres3];
    uint64_t* kmer_buffer1 = new uint64_t[(m_thres1) * span1];
    uint64_t* kmer_buffer2 = new uint64_t[(m_thres2) * span2];
    uint64_t* kmer_buffer3 = new uint64_t[(m_thres3) * span3];
    size_t no_minimiser1, no_minimiser2, no_minimiser3;
    uint64_t minimiser1, minimiser_rank1;
    uint64_t minimiser2, minimiser_rank2;
    uint64_t minimiser3, minimiser_rank3;
    size_t left_minimiser_position1, right_minimiser_position1;
    size_t left_minimiser_position2, right_minimiser_position2;
    size_t left_minimiser_position3, right_minimiser_position3;
    bool rolling1 = false;
    bool rolling2 = false;
    bool rolling3 = false;

    auto process = [&](auto&& rng)
    {
        for (auto&& kmer : rng) {
            if(!rolling1) {
                minimiser1 = find_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1);
                rolling1 = true;
            }
            else
                update_minimiser<1>(kmer.kmer_value, kmer.kmer_value_rev, minimiser1, left_minimiser_position1, right_minimiser_position1);

            if(minimiser1 == current_pos_minimiser1) {
                locate_buffer<1>(kmer_buffer1, offsets1, no_minimiser1, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1, positions, found_positions, found_kmers);
                rolling2 = false;
                rolling3 = false;
            }
            else if(minimiser1 != current_neg_minimiser1 && r1.contains(minimiser1, minimiser_rank1)) {
                const size_t minimiser_position1 = s1_select.select(minimiser_rank1);
                no_minimiser1 = s1_select.select(minimiser_rank1+1) - minimiser_position1;

                fill_buffer<1>(offsets1, kmer_buffer1, minimiser_position1, no_minimiser1, shift);
                locate_buffer<1>(kmer_buffer1, offsets1, no_minimiser1, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position1, right_minimiser_position1, positions, found_positions, found_kmers);
                current_pos_minimiser1 = minimiser1;
                rolling2 = false;
                rolling3 = false;
            }
            else {
                if(!rolling2) {
                    minimiser2 = find_minimiser<2>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2);
                    rolling2 = true;
                }
                else
                    update_minimiser<2>(kmer.kmer_value, kmer.kmer_value_rev, minimiser2, left_minimiser_position2, right_minimiser_position2);

                if (minimiser2 == current_pos_minimiser2) {
                    locate_buffer<2>(kmer_buffer2, offsets2, no_minimiser2, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2, positions, found_positions, found_kmers);
                    rolling3 = false;
                }
                else if(minimiser2 != current_neg_minimiser2 && r2.contains(minimiser2, minimiser_rank2)) {
                    const size_t minimiser_position2 = s2_select.select(minimiser_rank2);
                    no_minimiser2 = s2_select.select(minimiser_rank2+1) - minimiser_position2;

                    fill_buffer<2>(offsets2, kmer_buffer2, minimiser_position2, no_minimiser2, shift);
                    locate_buffer<2>(kmer_buffer2, offsets2, no_minimiser2, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position2, right_minimiser_position2, positions, found_positions, found_kmers);
                    current_pos_minimiser2 = minimiser2;
                    current_neg_minimiser1 = minimiser1;
                    rolling3 = false;
                }
                else {
                    if(!rolling3) {
                        minimiser3 = find_minimiser<3>(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position3, right_minimiser_position3);
                        rolling3 = true;
                    }
                    else
                        update_minimiser<3>(kmer.kmer_value, kmer.kmer_value_rev, minimiser3, left_minimiser_position3, right_minimiser_position3);

                    if(minimiser3 == current_pos_minimiser3) {
                        locate_buffer<3>(kmer_buffer3, offsets3, no_minimiser3, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position3, right_minimiser_position3, positions, found_positions, found_kmers);
                    }
                    else if(minimiser3 != current_neg_minimiser3 && r3.contains(minimiser3, minimiser_rank3)) {
                        const size_t minimiser_position3 = s3_select.select(minimiser_rank3);
                        no_minimiser3 = s3_select.select(minimiser_rank3+1) - minimiser_position3;

                        fill_buffer<3>(offsets3, kmer_buffer3, minimiser_position3, no_minimiser3, shift);
                        locate_buffer<3>(kmer_buffer3, offsets3, no_minimiser3, kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position3, right_minimiser_position3, positions, found_positions, found_kmers);
                        current_pos_minimiser3 = minimiser3;
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                    }
                    else {
                        if (auto it = hashmap.find(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev)); it != hashmap.end()) {
                            for(auto pos : it->second)
                                found_positions++;
                            found_kmers++;
                        }
                        current_neg_minimiser1 = minimiser1;
                        current_neg_minimiser2 = minimiser2;
                        current_neg_minimiser3 = minimiser3;
                    }
                }
            }
        }
    };

    if(shape != std::numeric_limits<uint32_t>::max())
        process(query | rshash::views::shapeview({.shape = shape}));
    else
        process(query | rshash::views::kmerview({.window_size = k}));

    delete[] kmer_buffer1;
    delete[] kmer_buffer2;
    delete[] kmer_buffer3;
    delete[] offsets1;
    delete[] offsets2;
    delete[] offsets3;

    return found_kmers;
}

#include "rshash.hpp"
#include "minimiser_views.hpp"
#include "shape_views.hpp"
#include "kxsort.h"

inline uint64_t mark_sequences(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t k,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &endpoints)
{
    size_t length = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        length += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << length << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "mark endpoints BV...\n";
    bit_vector sequences = bit_vector(length+33, 0);
    sequences[0] = 1;
    sequences[32] = 1;
    uint64_t j = 32;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), length+33);
    sequences = bit_vector();

    return kmers;
}

template<int level, typename MinimizerT>
inline void RSHash::filter_freq_minimizers(std::vector<MinimizerT> &minimizers,
    std::vector<uint8_t> &counts, size_t &no_minimizers, size_t &no_skmers)
{
    uint64_t threshold;
    if constexpr (level == 1)
        threshold = m_thres1;
    if constexpr (level == 2)
        threshold = m_thres2;
    if constexpr (level == 3)
        threshold = m_thres3;
    
    uint64_t current_minimizer = minimizers[0].first;
    uint64_t start = 0;
    size_t write_pos_idx = 0, write_min_idx = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i].first != current_minimizer) {
            size_t occurences = i - start;
            if(occurences <= threshold) {
                minimizers[write_min_idx++].first = current_minimizer;
                for (size_t j = start; j < i; j++)
                    minimizers[write_pos_idx++].second = minimizers[j].second;
                counts.push_back(static_cast<uint8_t>(occurences));
            }
            current_minimizer = minimizers[i].first;
            start = i;
        }
    }
    size_t occurences = minimizers.size() - start;
    if(minimizers.back().first == current_minimizer && occurences <= threshold) {
        minimizers[write_min_idx++].first = current_minimizer;
        for (size_t j = start; j < minimizers.size(); j++)
            minimizers[write_pos_idx++].second = minimizers[j].second;
        counts.push_back(static_cast<uint8_t>(occurences));
    }

    no_minimizers = write_min_idx;
    no_skmers = write_pos_idx;
}


template<int level, typename MinimizerT>
inline uint64_t RSHash::get_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers,
    std::vector<uint64_t> &unfreq_minimizers, std::vector<uint8_t> &counts)
{
    std::cout << std::bitset<64>(shape.value) << ' ' << window_size << " " << k << " " << overlap << " " << m1 << " " << m2 << " " << m3 << "\n";
    auto view = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    }();

    std::cout << "computing minimizers and positions...\n";
    std::vector<MinimizerT> minimizers;
    if constexpr (level == 1) {
        uint64_t length = 32;
        for(auto & sequence : sequences) {
            for(auto && minimizer : sequence | view)
                minimizers.emplace_back(MinimizerT{minimizer.minimiser_value, length + minimizer.range_position});
            length += sequence.size();
        }
    }
    else {
        for (auto skmerinfo : skmers) {
            auto skmer = sequences[skmerinfo.seq_id] | std::views::drop(skmerinfo.start) | std::views::take(skmerinfo.end - skmerinfo.start);
            size_t sequence_position = endpoints.select(skmerinfo.seq_id+1);    
            for (auto && minimizer : skmer | view)
                minimizers.emplace_back(MinimizerT{minimizer.minimiser_value, sequence_position + skmerinfo.start + minimizer.range_position});
        }
    }

    std::cout << "sorting minimizers...\n";
    // std::sort(minimizers.begin(), minimizers.end(), [](auto const& a, auto const& b) { return a.minimizer_value < b.minimizer_value; });
    // if(MinimizerT().second <= UINT32_MAX)
    //     kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer32());
    // else
        kx::radix_sort(minimizers.begin(), minimizers.end(), RadixTraitsMinimizer64());

    std::cout << "filtering frequent minimizers...\n";
    size_t no_minimizers, no_skmers;
    filter_freq_minimizers<level, MinimizerT>(minimizers, counts, no_minimizers, no_skmers);

    std::cout << "filling minimiser offsets...\n";
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, std::bit_width(endpoints.size()));
    for(size_t i = 0; i < no_skmers; i++)
        builder.push_back(minimizers[i].second);
    
    uint64_t* data_ptr = reinterpret_cast<uint64_t*>(minimizers.data());
    for(size_t i = 0; i < no_minimizers; i++)
        data_ptr[i] = minimizers[i].first;
    auto* raw = reinterpret_cast<uint64_t*>(minimizers.data());
    std::vector<uint64_t> tmp(raw, raw + no_minimizers);
    unfreq_minimizers.swap(tmp);

    minimizers.clear();

    if constexpr (level == 1)
        builder.build(offsets1);
    if constexpr (level == 2)  
        builder.build(offsets2);
    if constexpr (level == 3)
        builder.build(offsets3);

    return no_skmers;
}


inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding

    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : input | std::views::join | seqan3::views::to_rank) {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64) {
            packed.push_back(word);
            word = 0;
            shift = 0;
        }
    }

    if (shift != 0)
        packed.push_back(word);

    packed.push_back(UINT64_MAX); // padding

    return packed;
}


template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint8_t> &occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1;
    if constexpr (level == 2) return s2;
    if constexpr (level == 3) return s3;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 1)
        s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 2)
        s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 3)
        s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}

template<int level>
void RSHash::mark_occurences(const size_t total_occs, const std::vector<uint16_t> &occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 4) return s4;
    if constexpr (level == 5) return s5;
    }();

    s = bit_vector(total_occs+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < occurences.size(); i++) {
        j += occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 4)
        s4_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
    if constexpr (level == 5)
        s5_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), total_occs+1, 3);
}


template<int level>
size_t RSHash::get_frequent_skmers(
    const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const std::vector<SkmerInfo> &skmers, std::vector<SkmerInfo> &skmers_out)
{
    auto skmerview = [&]() {
        if constexpr (level == 1)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed1});
        if constexpr (level == 2)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed = seed2});
        if constexpr (level == 3)
            return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m3, .window_size = k, .seed = seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    if constexpr (level == 2) return r2;
    if constexpr (level == 3) return r3;
    }();
    
    size_t freq_kmers = 0;
    if constexpr (level == 1)
    {
        for(size_t sequence_id = 0; sequence_id < sequences.size(); ++sequence_id) {
            auto const & sequence = sequences[sequence_id];
            size_t start_position = 0;
            bool cur_freq, freq;
            uint64_t minimiser_rank;

            for(auto && minimiser : sequence | skmerview) {
                freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
                break;
            }
            for(auto && minimiser : sequence | skmerview) {
                cur_freq = r1.contains(minimiser.minimiser_value, minimiser_rank);
                if(freq && !cur_freq)
                    start_position = minimiser.range_position;
                if(!freq && cur_freq) {
                    skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, minimiser.range_position-1+k});
                    freq_kmers += minimiser.range_position - start_position;
                }
                freq = cur_freq;
            }
            if(!cur_freq) {
                skmers_out.emplace_back(SkmerInfo{sequence_id, start_position, sequence.size()});
                freq_kmers += sequence.size() - start_position - k + 1;
            }
        }
    }
    else
    {
        for(auto skmer : skmers) {
            auto skmer_view = sequences[skmer.seq_id] | std::views::drop(skmer.start) | std::views::take(skmer.end - skmer.start);
            bool cur_freq, freq;
            size_t start = 0;
            uint64_t minimiser_rank;

            for(auto && minimiser : skmer_view | skmerview) {
                freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                break;
            }
            for(auto && minimiser : skmer_view | skmerview) {
                cur_freq = r.contains(minimiser.minimiser_value, minimiser_rank);
                if(freq && !cur_freq)
                    start = minimiser.range_position;
                if(!freq && cur_freq) {
                    skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.start + minimiser.range_position-1+k});
                    freq_kmers += minimiser.range_position - start;
                }
                freq = cur_freq;
            }
            if(!cur_freq) {
                skmers_out.emplace_back(SkmerInfo{skmer.seq_id, skmer.start + start, skmer.end});
                freq_kmers += skmer.end - skmer.start - start - k + 1;
            }
        }
    }

    return freq_kmers;
}


void RSHash::fill_ht(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    const bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();

    if(loc) {
        gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts; // assert threshold <= 2^16
        if(threshold > 0) {
            for(const auto & skmer_info : freq_skmers) {
                size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
                size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
                auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
                for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                    if(use_shape) { // todo: symmetric shapes are canonical
                        uint64_t shape = _pext_u64(kmer.value, shape_mask);
                        uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                        if(kmer_counts[shape] < threshold)
                            kmer_counts[shape]++;
                        if(kmer_counts[shape_rc] < threshold)
                            kmer_counts[shape_rc]++;
                    }
                    else {
                        uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                        if(kmer_counts[canonical_kmer] < threshold)
                            kmer_counts[canonical_kmer]++;
                    }
                }
            }
        }
        FlatMapBuilder hashmap_builder, hashmap_rc_builder;
        for(const auto & skmer_info : freq_skmers) {
            size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
            size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
            auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
            const uint32_t skmer_pos = endpoints.select(skmer_info.seq_id+1) + skmer_info.start - overlap;
            uint32_t p = 0;
            for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                if(use_shape) { // todo: symmetric shapes are canonical
                    const uint64_t shape = _pext_u64(kmer.value, shape_mask);
                    const uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                    if(threshold > 0) {
                        if(kmer_counts[shape] < threshold)
                            hashmap_builder.add(shape, skmer_pos + p);
                        if(kmer_counts[shape_rc] < threshold)
                            hashmap_rc_builder.add(shape_rc, skmer_pos + p);
                    }
                    else {
                        hashmap_builder.add(shape, skmer_pos + p);
                        hashmap_rc_builder.add(shape_rc, skmer_pos + p);
                    }
                }
                else {
                    const uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                    if(threshold > 0) {
                        if(kmer_counts[canonical_kmer] < threshold)
                            hashmap_builder.add(canonical_kmer, skmer_pos + p);
                    }
                    else
                        hashmap_builder.add(canonical_kmer, skmer_pos + p); 
                }   
                p++;
            }
        }
        hashmap = std::move(hashmap_builder).build();
        hashmap_rc = std::move(hashmap_rc_builder).build();
    }
    else {
        for(const auto & skmer_info : freq_skmers) {
            size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
            size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
            auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
            // auto skmer = input[skmer_info.seq_id] | std::views::drop(skmer_info.start) | std::views::take(skmer_info.end - skmer_info.start);
            for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                if(use_shape) { // todo: symmetric shapes are canonical
                    const uint64_t shape = _pext_u64(kmer.value, shape_mask);
                    const uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                        hashset.insert(shape);
                        hashset_rc.insert(shape_rc);
                }
                else
                    hashset.insert(std::min<uint64_t>(kmer.value, kmer.value_rev));
            }
        }
    }
}


void RSHash::last_level(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input,
    const std::vector<SkmerInfo> &freq_skmers)
{
    const bool use_shape = shape.value != std::numeric_limits<uint32_t>::max();
    const uint64_t M = 1ULL << (2*window_size);

    if(loc) {
        gtl::flat_hash_map<uint64_t, uint16_t> kmer_counts; // assert threshold <= 2^16
        if(threshold > 0) {
            for(const auto & skmer_info : freq_skmers) {
                size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
                size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
                auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
                for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                    if(use_shape) { // todo: symmetric shapes are canonical
                        uint64_t shape = _pext_u64(kmer.value, shape_mask);
                        uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                        if(kmer_counts[shape] < threshold)
                            kmer_counts[shape]++;
                        if(kmer_counts[shape_rc] < threshold)
                            kmer_counts[shape_rc]++;
                    }
                    else {
                        uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                        if(kmer_counts[canonical_kmer] < threshold)
                            kmer_counts[canonical_kmer]++;
                    }
                }
            }
        }

        gtl::flat_hash_map<uint64_t, std::vector<uint64_t>> freq_kmers_map; // todo: 32 bit pos if possible
        gtl::flat_hash_map<uint64_t, std::vector<uint64_t>> freq_kmers_rc_map;
        for(const auto & skmer_info : freq_skmers) {
            size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
            size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
            auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
            const uint32_t skmer_pos = endpoints.select(skmer_info.seq_id+1) + skmer_info.start - overlap;
            // todo: rather save seq_id & skmer_pos in 64 bit
            uint32_t p = 0;
            for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                if(use_shape) { // todo: symmetric shapes are canonical
                    const uint64_t shape = _pext_u64(kmer.value, shape_mask);
                    const uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                    if(threshold > 0) {
                        if(kmer_counts[shape] < threshold)
                            freq_kmers_map[shape].push_back(skmer_pos + p);
                        if(kmer_counts[shape_rc] < threshold)
                            freq_kmers_rc_map[shape_rc].push_back(skmer_pos + p);
                    }
                    else {
                        freq_kmers_map[shape].push_back(skmer_pos + p);
                        freq_kmers_rc_map[shape_rc].push_back(skmer_pos + p);
                    }
                }
                else {
                    const uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                    if(threshold > 0) {
                        if(kmer_counts[canonical_kmer] < threshold)
                            freq_kmers_map[canonical_kmer].push_back(skmer_pos + p);
                    }
                    else
                        freq_kmers_map[canonical_kmer].push_back(skmer_pos + p);
                }   
                p++;
            }
        }
        kmer_counts.clear();

        std::vector<gtl::flat_hash_map<uint64_t, std::vector<uint64_t>>::const_iterator> order;
        order.reserve(freq_kmers_map.size());

        size_t total_positions = 0;
        for (auto it = freq_kmers_map.begin(); it != freq_kmers_map.end(); ++it) {
            order.push_back(it);
            total_positions += it->second.size();
        }

        std::sort(order.begin(), order.end(), [](auto a, auto b) { return a->first < b->first; });

        std::vector<uint64_t> freq_kmers;
        std::vector<uint16_t> occurrences;
        bits::compact_vector::builder builder;
        freq_kmers.reserve(order.size());
        occurrences.reserve(order.size());
        builder.resize(total_positions, std::bit_width(endpoints.size()));

        for (auto it : order) {
            freq_kmers.push_back(it->first);
            occurrences.push_back(static_cast<uint16_t>(it->second.size()));
            for (const auto& pos : it->second)
                builder.push_back(pos);
        }
        freq_kmers_map.clear();

        std::sort(freq_kmers.begin(), freq_kmers.end());
        r4 = sux::bits::EliasFano(freq_kmers, M);
        freq_kmers.clear();

        mark_occurences<4>(total_positions, occurrences);
        builder.build(offsets4);

        if(use_shape) { // todo: and shape not symmetric
            std::vector<gtl::flat_hash_map<uint64_t, std::vector<uint64_t>>::const_iterator> order2;
            order2.reserve(freq_kmers_rc_map.size());

            total_positions = 0;
            for (auto it = freq_kmers_rc_map.begin(); it != freq_kmers_rc_map.end(); ++it) {
                order2.push_back(it);
                total_positions += it->second.size();
            }

            std::sort(order2.begin(), order2.end(), [](auto a, auto b) { return a->first < b->first; });

            std::vector<uint64_t> freq_kmers_rc;
            std::vector<uint16_t> occurrences_rc;
            bits::compact_vector::builder builder_rc;
            freq_kmers_rc.reserve(order2.size());
            occurrences_rc.reserve(order2.size());
            builder_rc.resize(total_positions, std::bit_width(endpoints.size()));

            for (auto it : order2) {
                freq_kmers_rc.push_back(it->first);
                occurrences_rc.push_back(static_cast<uint16_t>(it->second.size()));
                for (const auto& pos : it->second)
                    builder_rc.push_back(pos);
            }
            freq_kmers_rc_map.clear();

            std::sort(freq_kmers_rc.begin(), freq_kmers_rc.end());
            r5 = sux::bits::EliasFano(freq_kmers_rc, M);
            freq_kmers_rc.clear();

            mark_occurences<5>(total_positions, occurrences_rc);
            builder_rc.build(offsets5);
        }
    }
    else {
        gtl::flat_hash_set<uint64_t> freq_kmers_map;
        gtl::flat_hash_set<uint64_t> freq_kmers_rc_map;
        for(const auto & skmer_info : freq_skmers) {
            size_t s = (skmer_info.start > overlap) ? skmer_info.start - overlap : 0;
            size_t e = std::min(skmer_info.end + overlap, input[skmer_info.seq_id].size()); // + overlap?
            auto skmer = input[skmer_info.seq_id] | std::views::drop(s) | std::views::take(e - s);
            for(auto && kmer : skmer | rshash::views::kmerview({.window_size = window_size})) {
                if(use_shape) {
                    const uint64_t shape = _pext_u64(kmer.value, shape_mask);
                    const uint64_t shape_rc = _pext_u64(kmer.value, shape_mask_rev);
                    if(shape != shape_rc) {
                        freq_kmers_map.insert(shape);
                        freq_kmers_rc_map.insert(shape_rc);
                    }
                    else
                        freq_kmers_map.insert(shape);
                }
                else {
                    const uint64_t canonical_kmer = std::min<uint64_t>(kmer.value, kmer.value_rev);
                    freq_kmers_map.insert(canonical_kmer);
                }
            }
        }
        std::vector<uint64_t> freq_kmers(freq_kmers_map.begin(), freq_kmers_map.end());
        freq_kmers_map.clear();
        std::sort(freq_kmers.begin(), freq_kmers.end());
        r4 = sux::bits::EliasFano(freq_kmers, M);
        freq_kmers.clear();

        std::vector<uint64_t> freq_kmers_rc(freq_kmers_rc_map.begin(), freq_kmers_rc_map.end());
        freq_kmers_rc_map.clear();
        std::sort(freq_kmers_rc.begin(), freq_kmers_rc.end());
        r5 = sux::bits::EliasFano(freq_kmers_rc, M);
        freq_kmers_rc.clear();
    }
}


void RSHash::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    const size_t text_length = endpoints.size();
    const size_t log_text_length = std::bit_width(text_length);

    std::vector<uint64_t> minimizers1;
    std::vector<uint8_t> minimizers1_occurences;
    std::vector<SkmerInfo> freq_skmers;
    uint64_t no_skmers1;
    if(log_text_length <= 32)
        no_skmers1 = get_minimizers<1, MinimizerInfo32>(input, freq_skmers, minimizers1, minimizers1_occurences);
    else
        no_skmers1 = get_minimizers<1, MinimizerInfo64>(input, freq_skmers, minimizers1, minimizers1_occurences);

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    mark_occurences<1>(no_skmers1, minimizers1_occurences);
    minimizers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    size_t freq_kmers;
    freq_kmers = get_frequent_skmers<1>(input, freq_skmers, freq_skmers);

    if(level > 1)
    {
        std::cout << "count minimizers2...\n";
        std::vector<uint64_t> minimizers2;
        std::vector<uint8_t> minimizers2_occurences;
        uint64_t no_skmers2;
        if(log_text_length <= 32)
            no_skmers2 = get_minimizers<2, MinimizerInfo32>(input, freq_skmers, minimizers2, minimizers2_occurences);
        else
            no_skmers2 = get_minimizers<2, MinimizerInfo64>(input, freq_skmers, minimizers2, minimizers2_occurences);

        std::cout << "build R_2...\n";
        const uint64_t M2 = 1ULL << (m2+m2);
        r2 = sux::bits::EliasFano(minimizers2, M2);
        minimizers2.clear();

        std::cout << "mark skmers2...\n";
        mark_occurences<2>(no_skmers2, minimizers2_occurences);
        minimizers2_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers2;
        freq_kmers = get_frequent_skmers<2>(input, freq_skmers, freq_skmers2);
        freq_skmers.clear();
        freq_skmers = freq_skmers2;
    }

    if(level > 2)
    {
        std::cout << "count minimizers3...\n";
        std::vector<uint64_t> minimizers3;
        std::vector<uint8_t> minimizers3_occurences;
        uint64_t no_skmers3;
        if(log_text_length <= 32)
            no_skmers3 = get_minimizers<3, MinimizerInfo32>(input, freq_skmers, minimizers3, minimizers3_occurences);
        else
            no_skmers3 = get_minimizers<3, MinimizerInfo64>(input, freq_skmers, minimizers3, minimizers3_occurences);

        std::cout << "build R_3...\n";
        const uint64_t M3 = 1ULL << (m3+m3);
        r3 = sux::bits::EliasFano(minimizers3, M3);
        minimizers3.clear();

        std::cout << "mark skmers3...\n";
        mark_occurences<3>(no_skmers3, minimizers3_occurences);
        minimizers3_occurences.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<SkmerInfo> freq_skmers3;
        freq_kmers = get_frequent_skmers<3>(input, freq_skmers, freq_skmers3);
        freq_skmers.clear();
        freq_skmers = freq_skmers3;
    }

    std::cout << "build final level...\n";
    if(use_ht)
        fill_ht(input, freq_skmers);
    else
        last_level(input, freq_skmers);

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}

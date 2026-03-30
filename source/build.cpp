#include "rshash.hpp"
#include "minimiser_views.hpp"


void RSHash::print_info() {
    const size_t N = text.size()*32;
    const size_t offset_width = std::bit_width(N);
    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);
    const uint64_t M3 = 1ULL << (m3+m3);
    const uint64_t no_minimizers1 = r1.rank(M1);
    const uint64_t no_skmers1 = s1.size();
    const uint64_t no_minimizers2 = r2.rank(M2);
    const uint64_t no_skmers2 = s2.size();
    const uint64_t no_minimizers3 = r3.rank(M3);
    const uint64_t no_skmers3 = s3.size();

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "textkmers: " << no_text_kmers <<  '\n';
    
    std::cout << "no minimiser1: " << no_minimizers1 << "\n";
    std::cout << "no distinct minimiser1: " << no_skmers1 << "\n";
    std::cout << "avg superkmers1: " << (double) no_skmers1/no_minimizers1 <<  '\n';
    std::cout << "no minimiser2: " << no_minimizers2 << "\n";
    std::cout << "no distinct minimiser2: " << no_skmers2 << "\n";
    std::cout << "avg superkmers2: " << (double) no_skmers2/no_minimizers2 <<  '\n';
    std::cout << "no minimiser3: " << no_minimizers3 << "\n";
    std::cout << "no distinct minimiser3: " << no_skmers3 << "\n";
    std::cout << "avg superkmers3: " << (double) no_skmers3/no_minimizers3 <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/no_text_kmers*100 << "%\n";

    std::cout << "density r1: " << (double) no_minimizers1/M1*100 << "%\n";
    std::cout << "density r2: " << (double) no_minimizers2/M2*100 << "%\n";
    std::cout << "density r3: " << (double) no_minimizers3/M3*100 << "%\n";
    std::cout << "density s1: " << (double) no_minimizers1/(no_skmers1+1)*100 <<  "%\n";
    std::cout << "density s2: " << (double) no_minimizers2/(no_skmers2+1)*100 <<  "%\n";
    std::cout << "density s3: " << (double) no_minimizers3/(no_skmers3+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/no_text_kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/no_text_kmers << "\n";
    std::cout << "offsets1: " << (double) no_skmers1*offset_width/no_text_kmers << "\n";
    std::cout << "offsets2: " << (double) no_skmers2*offset_width/no_text_kmers << "\n";
    std::cout << "offsets3: " << (double) no_skmers3*offset_width/no_text_kmers << "\n";
    std::cout << "Hashtable: " << (double) hashmap.capacity()*(sizeof(uint64_t) + 1)*8/no_text_kmers << "\n";
    std::cout << "R_1: " << (double) r1.bitCount()/no_text_kmers << "\n";
    std::cout << "R_2: " << (double) r2.bitCount()/no_text_kmers << "\n";
    std::cout << "R_3: " << (double) r3.bitCount()/no_text_kmers << "\n";
    std::cout << "S_1: " << (double) (no_skmers1+1)/no_text_kmers << "\n";
    std::cout << "S_2: " << (double) (no_skmers2+1)/no_text_kmers << "\n";
    std::cout << "S_3: " << (double) (no_skmers3+1)/no_text_kmers << "\n";
    
    std::cout << "total: " << (double) (no_skmers1*offset_width+no_skmers2*offset_width+no_skmers3*offset_width+2*N+r1.bitCount()+r2.bitCount()+r3.bitCount()+no_skmers1+1+s1_select.bitCount()+no_skmers2+1+s2_select.bitCount()+no_skmers3+1+s3_select.bitCount()+endpoints.bitCount()+65*hashmap.bucket_count())/no_text_kmers << "\n";
}


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


inline uint64_t get_unfrequent_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const size_t minimiser_length, const size_t threshold, const size_t k, const uint64_t seed,
    std::vector<uint64_t> &unfreq_minimizers, std::vector<uint8_t> &counts)
{
    std::vector<uint64_t> minimizers;

    std::cout << "computing minimizers...\n";
    auto minimiserview = rshash::views::xor_minimiser_and_positions({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimiser : sequence | minimiserview)
            minimizers.emplace_back(minimiser.minimiser_value);
    }

    std::cout << "sorting minimizers...\n";
    std::sort(minimizers.begin(), minimizers.end());

    std::cout << "get unfrequent minimizers...\n";
    uint64_t current_minimizer = minimizers[0];
    uint64_t occurences = 1;
    uint64_t no_skmers = 0;
    for(size_t i = 1; i < minimizers.size(); i++) {
        if(minimizers[i] != current_minimizer) {
            if(occurences <= threshold) {
                unfreq_minimizers.emplace_back(current_minimizer);
                counts.emplace_back(static_cast<uint8_t>(occurences));
                no_skmers += occurences;
            }
            current_minimizer = minimizers[i];
            occurences = 1;
        }
        else
            occurences++;
    }
    if(minimizers.back() == current_minimizer && occurences <= threshold) {
        unfreq_minimizers.emplace_back(current_minimizer);
        counts.emplace_back(static_cast<uint8_t>(occurences));
        no_skmers += occurences;
    }

    minimizers.clear();

    return no_skmers;
}



inline std::vector<uint64_t> pack_dna4_to_uint64(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    auto ranks = input | std::views::join | seqan3::views::to_rank;

    std::vector<uint64_t> packed;
    packed.push_back(UINT64_MAX); // padding

    uint64_t word = 0;
    size_t shift = 0;

    for (uint8_t r : ranks)
    {
        word |= uint64_t(r) << shift; // pack 2 bits per base
        shift += 2;

        if (shift == 64)
        {
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
void RSHash::mark_minimizer_occurences(const size_t no_skmers, const std::vector<uint8_t> &minimizer_occurences)
{
    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1;
    else if constexpr (level == 2) return s2;
    else if constexpr (level == 3) return s3;
    }();

    s = bit_vector(no_skmers+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < minimizer_occurences.size(); i++) {
        j += minimizer_occurences[i];
        s[j] = 1;
    }

    if constexpr (level == 1)
        s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);
    if constexpr (level == 2)
        s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);
    if constexpr (level == 3)
        s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);
}


template<int level>
void RSHash::fill_minimizer_offsets(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    std::vector<size_t> &skmer_positions, std::vector<uint8_t> &minimizer_occurences,
    const size_t text_length, const size_t no_skmers)
{
    auto view = [&]() {
    if constexpr (level == 1)
        return rshash::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    else if constexpr (level == 2)
        return rshash::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    else if constexpr (level == 3)
        return rshash::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    else if constexpr (level == 2) return r2;
    else if constexpr (level == 3) return r3;
    }();

    auto & s = [&]() -> auto& {
    if constexpr (level == 1) return s1_select;
    else if constexpr (level == 2) return s2_select;
    else if constexpr (level == 3) return s3_select;
    }();

    const size_t offset_width = std::bit_width(text_length);
    bits::compact_vector::builder builder;
    builder.resize(no_skmers, offset_width);

    std::fill(minimizer_occurences.begin(), minimizer_occurences.end(), 0);

    size_t length = 32;
    size_t skmer_idx = 0;
    for(auto & sequence : sequences) {
        for (auto && minimiser : sequence | view) {
            if(uint64_t minimizer_rank = r.rank(minimiser.minimiser_value); r.rank(minimiser.minimiser_value+1)-minimizer_rank) {
                size_t minimizer_idx = s.select(minimizer_rank);
                if constexpr (level == 1)
                    builder.set(minimizer_idx + minimizer_occurences[minimizer_rank], length + minimiser.range_position);
                else
                    builder.set(minimizer_idx + minimizer_occurences[minimizer_rank], skmer_positions[skmer_idx] + minimiser.range_position);
                minimizer_occurences[minimizer_rank]++;
            }
        }
        skmer_idx++;
        length += sequence.size();
    }

    if constexpr (level == 1)
        builder.build(offsets1);
    else if constexpr (level == 2)
        builder.build(offsets2);
    else if constexpr (level == 3)
        builder.build(offsets3);
}


template<int level>
void RSHash::get_frequent_skmers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences, std::vector<size_t> &positions,
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &freq_skmers, std::vector<size_t> &skmer_positions)
{
    auto skmerview = [&]() {
    if constexpr (level == 1)
        return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed = seed1});
    else if constexpr (level == 2)
        return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed = seed2});
    else if constexpr (level == 3)
        return rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m3, .window_size = k, .seed = seed3});
    }();

    auto & r = [&]() -> auto& {
    if constexpr (level == 1) return r1;
    else if constexpr (level == 2) return r2;
    else if constexpr (level == 3) return r3;
    }();
    
    size_t length = 32;
    size_t skmer_idx = 0;
    for(auto & sequence : sequences) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | skmerview) {
            freq = r.rank(minimiser.minimiser_value+1)-r.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            cur_freq = r.rank(minimiser.minimiser_value+1)-r.rank(minimiser.minimiser_value);
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                seqan3::bitpacked_sequence<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers.emplace_back(skmer);
                if constexpr (level == 1)
                    skmer_positions.emplace_back(length + start_position);
                else
                    skmer_positions.emplace_back(positions[skmer_idx] + start_position);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            seqan3::bitpacked_sequence<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers.emplace_back(skmer);
            if constexpr (level == 1)
                skmer_positions.emplace_back(length + start_position);
            else
                skmer_positions.emplace_back(positions[skmer_idx] + start_position);
        }
        length += sequence.size();
        skmer_idx++;
    }
}


void RSHash::build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>& input)
{
    no_text_kmers = mark_sequences(input, k, endpoints);
    size_t text_length = endpoints.size();

    std::vector<uint64_t> minimizers1;
    std::vector<uint8_t> minimizers1_occurences;
    const uint64_t no_skmers1 = get_unfrequent_minimizers(input, m1, m_thres1, k, seed1, minimizers1, minimizers1_occurences);
    const size_t no_minimizers1 = minimizers1.size();

    std::cout << "build R_1...\n";
    const uint64_t M1 = 1ULL << (m1+m1);
    r1 = sux::bits::EliasFano(minimizers1, M1);
    minimizers1.clear();

    std::cout << "mark skmers1...\n";
    mark_minimizer_occurences<1>(no_skmers1, minimizers1_occurences);

    std::cout << "filling offsets_1...\n";
    std::vector<size_t> skmer_positions;
    fill_minimizer_offsets<1>(input, skmer_positions, minimizers1_occurences, text_length, no_skmers1);
    minimizers1_occurences.clear();

    std::cout << "get frequent skmers...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers;
    get_frequent_skmers<1>(input, skmer_positions, freq_skmers, skmer_positions);

    if(level > 1)
    {
        std::cout << "count minimizers2...\n";
        std::vector<uint64_t> minimizers2;
        std::vector<uint8_t> minimizers2_occurences;
        const uint64_t no_skmers2 = get_unfrequent_minimizers(freq_skmers, m2, m_thres2, k, seed2, minimizers2, minimizers2_occurences);
        const size_t no_minimizers2 = minimizers2.size();

        std::cout << "build R_2...\n";
        std::sort(minimizers2.begin(), minimizers2.end());
        const uint64_t M2 = 1ULL << (m2+m2);
        r2 = sux::bits::EliasFano(minimizers2, M2);

        std::cout << "mark skmers2...\n";
        mark_minimizer_occurences<2>(no_skmers2, minimizers2_occurences);

        std::cout << "filling offsets_2...\n";
        fill_minimizer_offsets<2>(freq_skmers, skmer_positions, minimizers2_occurences, text_length, no_skmers2);
        minimizers2_occurences.clear();
        skmer_positions.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers_;
        std::vector<size_t> skmer_positions_;
        get_frequent_skmers<2>(freq_skmers, skmer_positions, freq_skmers_, skmer_positions_);
        freq_skmers = freq_skmers_;
        skmer_positions = skmer_positions_;
    }

    if(level > 2)
    {
        std::cout << "count minimizers3...\n";
        std::vector<uint64_t> minimizers3;
        std::vector<uint8_t> minimizers3_occurences;
        const uint64_t no_skmers3 = get_unfrequent_minimizers(freq_skmers, m3, m_thres3, k, seed3, minimizers3, minimizers3_occurences);
        const size_t no_minimizers3 = minimizers3.size();

        std::cout << "build R_3...\n";
        std::sort(minimizers3.begin(), minimizers3.end());
        const uint64_t M3 = 1ULL << (m3+m3);
        r3 = sux::bits::EliasFano(minimizers3, M3);

        std::cout << "mark skmers3...\n";
        mark_minimizer_occurences<3>(no_skmers3, minimizers3_occurences);

        std::cout << "filling offsets_3...\n";
        fill_minimizer_offsets<3>(freq_skmers, skmer_positions, minimizers3_occurences, text_length, no_skmers3);
        minimizers3_occurences.clear();
        skmer_positions.clear();

        std::cout << "get frequent skmers...\n";
        std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> freq_skmers_;
        std::vector<size_t> skmer_positions_;
        get_frequent_skmers<3>(freq_skmers, skmer_positions, freq_skmers_, skmer_positions_);
        freq_skmers = freq_skmers_;
        skmer_positions = skmer_positions_;
    }

    std::cout << "build HT...\n";
    for(auto & sequence : freq_skmers)
        for(auto && kmer : sequence | rshash::views::kmerview({.window_size = k}))
            hashmap.insert(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev));

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    print_info();
}


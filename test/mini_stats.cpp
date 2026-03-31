#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <sharg/all.hpp>
#include <sux/bits/SimpleSelect.hpp>

#include "../source/EliasFano.hpp"
#include "../source/minimiser_views.hpp"

using namespace seqan3::contrib::sdsl;

struct cmd_arguments {
    std::filesystem::path i{};
    uint8_t k{31};
    uint8_t m{18};
    uint8_t t{64};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "kmer", .description = "k-mer length"});
    parser.add_option(args.m, sharg::config{.short_id = 'm', .long_id = "mini", .description = "minimiser length"});
    parser.add_option(args.t, sharg::config{.short_id = 't', .long_id = "thres", .description = "threshold"});
}


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}


uint64_t count_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const uint8_t minimiser_length, const uint8_t threshold, const uint8_t k, const uint64_t seed,
    std::vector<uint64_t> &unfreq_minimizers, std::vector<uint8_t> &counts, uint64_t &kmers)
{
    // std::cout << "computing minimizers...\n";
    std::vector<uint64_t> minimizers;
    auto minimiserview = rshash::views::xor_minimiser_and_positions({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimiser : sequence | minimiserview)
            minimizers.emplace_back(minimiser.minimiser_value);
        kmers += sequence.size() - k + 1;
    }

    // std::cout << "sorting minimizers...\n";
    std::sort(minimizers.begin(), minimizers.end());

    // std::cout << "counting minimizers...\n";
    uint64_t current_minimizer = minimizers[0];
    uint64_t occurences = 0;
    uint64_t no_skmers = 0;
    for(uint64_t minimizer : minimizers) {
        if(minimizer != current_minimizer) {
            if(occurences <= threshold) {
                unfreq_minimizers.emplace_back(current_minimizer);
                counts.emplace_back(static_cast<uint8_t>(occurences));
                no_skmers += occurences;
            }
            current_minimizer = minimizer;
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


std::vector<uint64_t> get_uncovered_kmers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const size_t m, const size_t t, const size_t k, const uint64_t seed,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> &r, sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s)
{
    auto skmerview = rshash::views::xor_minimiser_and_skmer_positions({.minimiser_size = m, .window_size = k, .seed = seed});
    
    std::vector<uint64_t> kmers_uncovered(t, 0);
    for(auto & sequence : sequences) {
        std::vector<bool> cur_mini_unfreq(t), mini_unfreq(t, false);
        std::vector<size_t> start_positions(t, 0);
        uint64_t minimiser_rank;

        for(auto && minimiser : sequence | skmerview) {
            if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                for(size_t i=0; i < t; i++)
                    mini_unfreq[i] = frequency <= i+1;
            }
            break;
        }
        for(auto && minimiser : sequence | skmerview) {
            if(minimiser_rank = r.rank(minimiser.minimiser_value), r.rank(minimiser.minimiser_value + 1) - minimiser_rank) {
                size_t frequency = s.select(minimiser_rank+1)-s.select(minimiser_rank);
                
                for(size_t i=0; i < t; i++) {
                    cur_mini_unfreq[i] = frequency <= i+1;
                    if(mini_unfreq[i] && !cur_mini_unfreq[i])
                        start_positions[i] = minimiser.range_position;
                    if(!mini_unfreq[i] && cur_mini_unfreq[i])
                        kmers_uncovered[i] += minimiser.range_position - start_positions[i];
                    mini_unfreq[i] = cur_mini_unfreq[i];
                }
            }
            else {
                for(size_t i=0; i < t; i++) {
                    if(mini_unfreq[i])
                        start_positions[i] = minimiser.range_position;
                    mini_unfreq[i] = false;
                }
            }
        }
        for(size_t i=0; i < t; i++) {
            if(!cur_mini_unfreq[i])
                kmers_uncovered[i] += sequence.size() - k + 1 - start_positions[i];
        }
        
    }

    return kmers_uncovered;
}


void stats(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t m, const size_t t, const size_t k)
{
    uint64_t kmers = 0;
    std::vector<uint64_t> minimizers;
    std::vector<uint8_t> counts;
    uint64_t no_skmers = count_minimizers(input, m, t, k, 1, minimizers, counts, kmers);
    uint64_t no_minimizers = minimizers.size();

    const uint64_t M = 1ULL << (m+m);
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r = sux::bits::EliasFano(minimizers, M);
    minimizers.clear();

    bit_vector s = bit_vector(no_skmers+1, 0);
    s[0] = 1;
    uint64_t j = 0;
    for(size_t i = 0; i < counts.size(); i++) {
        j += counts[i];
        s[j] = 1;
    }
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), no_skmers+1, 3);

    // std::cout << "compute kmer coverage...\n";
    std::vector<uint64_t> kmers_uncovered = get_uncovered_kmers(input, m, t, k, 1, r, s_select);

    uint64_t counter[t] = {0};
    for(uint8_t count : counts) {
        counter[count-1]++;
    }

    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    std::cout << "frequency, superkmers, kmers covered\n";
    for(uint8_t j=0; j < t; j++) {
        cum += counter[j];
        cum_skmers += (j+1)*counter[j];
        std::cout << (j+1) << ',' << std::setprecision(8) << (double) counter[j]/no_minimizers*100 << ',' << std::setprecision(8) << (double) (kmers-kmers_uncovered[j])/kmers*100 << '\n';
    }
    
}


int main(int argc, char** argv)
{
    sharg::parser parser{"ministats", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);
    stats(text, args.m, args.t, args.k);

    return 0;
}
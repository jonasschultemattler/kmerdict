#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <sharg/all.hpp>


#include "../source/minimiser_views.hpp"


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


void count_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const size_t minimiser_length, const size_t threshold, const size_t k, const uint64_t seed,
    std::vector<uint8_t> &counts)
{
    std::vector<uint64_t> minimizers;

    // std::cout << "computing minimizers...\n";
    auto minimiserview = rshash::views::xor_minimiser_and_positions({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimiser : sequence | minimiserview)
            minimizers.emplace_back(minimiser.minimiser_value);
    }

    // std::cout << "sorting minimizers...\n";
    std::sort(minimizers.begin(), minimizers.end());

    // std::cout << "counting minimizers...\n";
    uint64_t current_minimizer = minimizers[0];
    uint64_t occurences = 1;
    for(uint64_t minimizer : minimizers) {
        if(minimizer != current_minimizer) {
            if(occurences <= threshold)
                counts.emplace_back(static_cast<uint8_t>(occurences));
            current_minimizer = minimizer;
            occurences = 1;
        }
        else
            occurences++;
    }
    if(minimizers.back() == current_minimizer && occurences <= threshold)
        counts.emplace_back(static_cast<uint8_t>(occurences));

    minimizers.clear();
}


void stats(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input, const size_t m, const size_t t, const size_t k)
{
    std::vector<uint8_t> counts;
    count_minimizers(input, m, t, k, 1, counts);

    uint64_t minimizers = counts.size();
    uint64_t skmers = 0;
    uint64_t counter[t] = {0};
    for(uint8_t count : counts) {
        counter[count-1]++;
        skmers += count;
    }

    uint64_t cum = 0;
    uint64_t cum_skmers = 0;
    for(uint8_t j=0; j < t; j++) {
        cum += counter[j];
        cum_skmers += (j+1)*counter[j];
        std::cout << (j+1) << ',' << counter[j] << '\n';
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
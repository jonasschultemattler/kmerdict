#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

using namespace seqan3::contrib::sdsl;

struct cmd_arguments {
    std::filesystem::path i{};
    std::filesystem::path o{};
    uint8_t k{31};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.o, sharg::config{.short_id = 'o', .long_id = "output", .description = "provide output file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
    if(!parser.is_option_set('o'))
        throw sharg::user_input_error("provide output file.");
    if(!parser.is_option_set('k'))
        throw sharg::user_input_error("specify k");
    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna5;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna5>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna5> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}


inline uint64_t preprocess(
    std::vector<seqan3::bitpacked_sequence<seqan3::dna5>> const & input,
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> & output,
    const uint8_t k)
{
    uint64_t valid_length = 0;
    uint64_t valid_kmers = 0;
    uint64_t invalid_kmers = 0;
    uint64_t N = 0;

    for (auto const & record : input) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        uint64_t valid_run = 0;
        for (auto const nuc : record) {
            if (seqan3::to_rank(nuc) == 4) {
                ++invalid_kmers;

                if (valid_run >= k) {
                    valid_kmers += valid_run - k + 1;
                    output.push_back(std::move(seq));
                }

                valid_run = 0;
                seq = seqan3::bitpacked_sequence<seqan3::dna4>{};
            }
            else {
                seq.push_back(seqan3::dna4{}.assign_rank(seqan3::to_rank(nuc)));
                ++valid_length;
                ++valid_run;
            }
        }
        if (valid_run >= k) {
            valid_kmers += valid_run - k + 1;
            output.push_back(std::move(seq));
        }
        N += record.size();
    }

    std::cout << "no sequences: " << input.size() << "\n";
    std::cout << "input text length: " << N << "\n";
    std::cout << "valid text length: " << valid_length << "\n";
    std::cout << "valid kmers: " << valid_kmers << "\n";
    std::cout << "invalid kmers: " << invalid_kmers << " " << (double) invalid_kmers / valid_length * 100 << "%\n";
    std::cout << "output sequences: " << output.size() << "\n";

    return N;
}

int main(int argc, char** argv)
{
    sharg::parser parser{"parser", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna5>> input;
    load_file(args.i, input);

    std::cout << "parsing sequences...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> output;
    preprocess(input, output, args.k);

    std::cout << "saving output...\n";
    {
        auto stream = seqan3::sequence_file_output{args.o};
        for (size_t i = 0; i < output.size(); ++i) {
            stream.emplace_back(output[i], "seq_" + std::to_string(i));
        }
    }
 
    return 0;
}
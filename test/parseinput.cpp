#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include "../source/util.hpp"
#include "../source/EliasFano.hpp"
#include <gtl/phmap.hpp>

using namespace seqan3::contrib::sdsl;

struct cmd_arguments {
    std::filesystem::path i{};
    uint8_t k{31};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
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


inline uint64_t parse_sequences(
    std::vector<seqan3::bitpacked_sequence<seqan3::dna5>> const & input,
    size_t const k,
    std::vector<uint64_t> & packed,
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> & endpoints)
{
    uint64_t valid_length = 0;
    uint64_t valid_kmers = 0;
    uint64_t invalid_kmers = 0;
    uint64_t N = 0;

    for (auto const & record : input) {
        uint64_t valid_run = 0;
        for (auto const nuc : record) {
            if (seqan3::to_rank(nuc) == 4) {
                ++invalid_kmers;

                if (valid_run >= k)
                    valid_kmers += valid_run - k + 1;

                valid_run = 0;
            }
            else {
                ++valid_length;
                ++valid_run;
            }
        }
        if (valid_run >= k)
            valid_kmers += valid_run - k + 1;
        N += record.size();
    }

    std::cout << "no sequences: " << input.size() << "\n";
    std::cout << "valid text length: " << valid_length << "\n";
    std::cout << "valid kmers: " << valid_kmers << "\n";
    std::cout << "invalid kmers: " << invalid_kmers << " " << (double) invalid_kmers / valid_length * 100 << "%\n";
    

    packed.push_back(UINT64_MAX);
    
    std::vector<size_t> sequences;
    sequences.push_back(0);
    sequences.push_back(32);

    uint64_t pos = 32;
    uint64_t word = 0;
    size_t shift = 0;

    for (auto const & record : input)
    {
        bool valid = true;
        for (uint8_t const rank : record | seqan3::views::to_rank) {
            if(rank == 4) {
                valid = false;
                continue;
            }

            if(!valid) {
                sequences.push_back(pos);
                valid = true;
            }

            word |= (uint64_t{rank} << shift);
            shift += 2;

            ++pos;

            if(shift == 64) {
                packed.push_back(word);
                word = 0;
                shift = 0;
            }
        }

        sequences.push_back(pos);
    }

    if(shift != 0)
        packed.push_back(word);

    packed.push_back(UINT64_MAX);

    endpoints = sux::bits::EliasFano(sequences, pos + 33);
    sequences.clear();

    std::cout << "text BV length/2: " << packed.size()*32 << "\n";
    std::cout << "endpoints pos: " << pos - 32 << "\n";
    std::cout << "no valid sequences: " << endpoints.rank(pos+1) << "\n";
    std::cout << "text length: " << endpoints.select(endpoints.rank(pos+1)) << "\n";

    return N;
}


const inline uint64_t get_word64(std::vector<uint64_t> const & text, uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t get_base(std::vector<uint64_t> const & text, uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
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
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints = sux::bits::EliasFano<sux::util::AllocType::MALLOC>(std::vector<uint64_t>{}, 1);
    std::vector<uint64_t> text;
    parse_sequences(input, args.k, text, endpoints);

    std::cout << "text vector length: " << text.size()*32 << "\n";
    size_t text_length = endpoints.select(endpoints.rank((text.size()-1)*32)) - 32;
    
    std::cout << "text length: " << text_length << "\n";
    std::cout << "text length: " << endpoints.select(1444+1) << "\n";

    const uint64_t kmer_mask = compute_mask(args.k);
    const uint64_t shift = (args.k-1)*2;

    gtl::flat_hash_map<uint64_t, uint32_t> ht;

    uint64_t kmer = get_word64(text, 32) & kmer_mask;
    uint64_t kmer_rc = crc(kmer, args.k);
    for(size_t i = 32+args.k; i < (text.size()-1)*32; i++) {
        // uint64_t kmer = get_word64(text, i) & kmer_mask;
        // uint64_t kmer_rc = crc(kmer, args.k);
        uint64_t base = get_base(text, i);
        // kmer = (kmer >> 2) | (base << shift);
        kmer = ((kmer << 2) | base) & kmer_mask;
        kmer_rc = kmer_rc >> 2 | ((base ^ 0b11) << shift);
        // kmer_rc = ((kmer_rc << 2) | (base ^ 0b11)) & kmer_mask;
        uint64_t canonical_kmer = std::min<uint64_t>(kmer, kmer_rc);
        // ht[kmer_rc]++;
        ht[kmer]++;
    }

    std::cout << "number of distinct kmers: " << ht.size() << "\n";
 
    return 0;
}
#include "include/benchmarks.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void usage(const char *exe)
{
    std::cerr << "Usage:\n"
              << "  " << exe << " --matrix PATH --output DIR [--repeats 5] [--spy]\n"
              << "  " << exe << " --thomas-scale N --output DIR [--repeats 5]\n"
              << "  " << exe << " --dense-scale N --output DIR [--repeats 5]\n";
}

std::size_t number(const std::string& str)
{
    std::size_t              end   = 0;
    const unsigned long long value = std::stoull(str, &end);
    if (end != str.size() || value == 0)
        throw std::invalid_argument("Expected a positive integer: " + str);
    return static_cast<std::size_t>(value);
}
}  // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc == 1)
        {
            usage(argv[0]);
            return EXIT_FAILURE;
        }

        std::string           mode, matrix_path;
        std::size_t           n = 0;
        std::filesystem::path output;
        benchmarks::Options   options;

        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--spy")
            {
                options.spy_plots = true;
                continue;
            }
            if (arg == "--help" || arg == "-h")
            {
                usage(argv[0]);
                return EXIT_SUCCESS;
            }
            if (i + 1 >= argc)
                throw std::invalid_argument("Missing value for " + arg);
            const std::string value = argv[++i];
            if (arg == "--matrix")
            {
                mode        = "matrix";
                matrix_path = value;
            }
            else if (arg == "--thomas-scale")
            {
                mode = "thomas";
                n    = number(value);
            }
            else if (arg == "--dense-scale")
            {
                mode = "dense";
                n    = number(value);
            }
            else if (arg == "--output")
                output = value;
            else if (arg == "--repeats")
                options.repeats = number(value);
            else
                throw std::invalid_argument("Unknown option: " + arg);
        }
        if (output.empty())
            throw std::invalid_argument("--output is required");
        if (options.repeats < 5)
            throw std::invalid_argument("At least five repeats required");
        if (mode == "matrix" && !matrix_path.empty())
            benchmarks::run_matrix(matrix_path, output, options);
        else if (mode == "thomas" && n > 0)
            benchmarks::run_thomas_scale(n, output, options);
        else if (mode == "dense" && n > 0)
            benchmarks::run_dense_scale(n, output, options);
        else
            throw std::invalid_argument("Select a valid mode");
        return EXIT_SUCCESS;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}

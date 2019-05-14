#pragma once
#include <iostream>
#include <stdio.h>
#include "cli11/CLI11.hpp"
#include "parameters.h"

// Failure message function for CLI parser
inline std::string failure_message(const CLI::App* parser, const CLI::Error& error) {

    std::string message = "";

    if (error.what() == std::string("A subcommand is required")) {
        message = "\nSubcommand error: missing or invalid subcommand\n\n" + parser->help();
    } else if (error.get_exit_code() == 106) {  // 106 corresponds to wrong argument type
        message = "\nArgument error: " + std::string(error.what()) + "\n\n" + parser->help();
    } else {
        message = "\nError: " + std::string(error.what()) + "\n\n" + parser->help();
    }

    return message;
}


// Formatter for CLI
class CustomFormatter : public CLI::Formatter {

    public:

        uint column_widths[3] {0, 0, 0};  // Will be used to store the maximum width of each column : flags, type, description
        uint border_width = 4;  // Space between two columns

        // Formatter for an Option line, overrides the same function from CLI::Formatter
        virtual std::string make_option(const CLI::Option* opt, bool is_positional) const {

            std::string option = "", name = "", type = "", description = "", default_value = "", required = "REQUIRED";

            // Generate option name, if positional -> just the name, if not positional -> <short_flag, long_flag>
            is_positional ? name = opt->get_name()[0] : name = "-" + opt->get_snames()[0] + ", --" + opt->get_lnames()[0];
            type = opt->get_type_name();
            description = opt->get_description();
            default_value = opt->get_defaultval();

            // Generate the help string for this option, adding the right number of spaces after each column based on column_widths
            option = name + std::string(border_width + column_widths[0] - name.size(), ' ');
            option += type + std::string(border_width + column_widths[1] - type.size(), ' ');
            option += description + std::string(border_width + column_widths[2] - description.size(), ' ');
            if (opt->get_required()) default_value = required;
            if (default_value != "") option += "[" + default_value + "]";
            option += "\n";

            return option;
        }

        void set_column_widths(CLI::App& parser) {
            std::string tmp = "";
            for (auto opt: parser.get_options()) {
                opt->get_positional() ? tmp = opt->get_name()[0] : tmp = "-" + opt->get_snames()[0] + ", --" + opt->get_lnames()[0];
                if (tmp.size() > this->column_widths[0]) this->column_widths[0] = static_cast<uint>(tmp.size());
                tmp = opt->get_type_name();
                if (tmp.size() > this->column_widths[1]) this->column_widths[1] = static_cast<uint>(tmp.size());
                tmp = opt->get_description();
                if (tmp.size() > this->column_widths[2]) this->column_widths[2] = static_cast<uint>(tmp.size());
            }
        }

};


// Argument parsing main function
inline Parameters parse_args(int& argc, char** argv) {

    CLI::App parser {""};  // Parser instance from CLI App parser
    Parameters parameters;

    std::shared_ptr<CustomFormatter> formatter(new CustomFormatter);

    // Main parser options
    parser.formatter(formatter);  // Set custom help format defined above
    parser.require_subcommand();  // Check that there is a subcommand
    parser.failure_message(failure_message);  // Formatting for error message

    CLI::App* stats = parser.add_subcommand("stats", "Compute stats from MAF file");

    CLI::Option* option = stats->add_option("-m, --maf-file", parameters.maf_file_path, "Path to a MAF file");
    option->required();
    option->check(CLI::ExistingFile);

    stats->add_option("-a, --alignability-wig", parameters.alignability_wig_file_path, "Path to an wig output file for alignability results");
    stats->add_option("-A, --alignability-table", parameters.alignability_table_file_path, "Path to an output file for alignability results");
    stats->add_option("-i, --identity-wig", parameters.identity_wig_file_path, "Path to an wig output file for identity results");
    stats->add_option("-l, --major-allele-wig", parameters.major_allele_wig_file_path, "Path to an wig output file for major allele results");
    stats->add_option("-c, --allele-count-wig", parameters.allele_count_wig_file_path, "Path to an wig output file for allele count results");

    // The parser throws an exception upon failure and implements an exit() method which output an error message and returns the exit code.
    try {

        parser.parse(argc, argv);

    } catch (const CLI::ParseError &e) {

        if (parser.get_subcommands().size() > 0) {

            std::string tmp = "";
            for (auto opt: parser.get_subcommands()[0]->get_options()) {
                opt->get_positional() ? tmp = opt->get_name()[0] : tmp = "-" + opt->get_snames()[0] + ", --" + opt->get_lnames()[0];
                if (tmp.size() > formatter->column_widths[0]) formatter->column_widths[0] = static_cast<uint>(tmp.size());
                tmp = opt->get_type_name();
                if (tmp.size() > formatter->column_widths[1]) formatter->column_widths[1] = static_cast<uint>(tmp.size());
                tmp = opt->get_description();
                if (tmp.size() > formatter->column_widths[2]) formatter->column_widths[2] = static_cast<uint>(tmp.size());
            }

        } else {
            formatter->column_widths[0] = 10;
            formatter->column_widths[1] = 0;
            formatter->column_widths[2] = 50;
        }

        exit(parser.exit(e));

    }

    CLI::App* subcommand = parser.get_subcommands()[0];
    parameters.command = subcommand->get_name();

    if (subcommand->count("-a") + subcommand->count("-A") + subcommand->count("-i") + subcommand->count("-l") + subcommand->count("-c") == 0) {
        formatter->set_column_widths(*subcommand);
        std::cerr << "\nArgument error : at least one output file is required\n\n";
        std::cerr << parser.help();
        exit(1);
    }

    return parameters;
}

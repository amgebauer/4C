// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_dat_file_utils.hpp"

#include "4C_io_inputreader.hpp"

FOUR_C_NAMESPACE_OPEN

void Core::IO::DatFileUtils::print_section_header(std::ostream& out, const std::string& header)
{
  constexpr std::size_t max_line_width = 65ul;
  FOUR_C_THROW_UNLESS(header.length() <= max_line_width, "Header '%s' too long", header.c_str());

  out << "--";
  out << std::string(std::max(max_line_width - header.length(), 0ul), '-');
  out << header << '\n';
}



void Core::IO::DatFileUtils::print_section(std::ostream& out, const std::string& header,
    const std::vector<Input::LineDefinition>& possible_lines)
{
  print_section_header(out, header);

  for (const auto& line : possible_lines)
  {
    out << "// ";
    line.print(out);
    out << '\n';
  }
}


std::vector<Input::LineDefinition> Core::IO::DatFileUtils::read_all_lines_in_section(
    Core::IO::DatFileReader& reader, const std::string& section,
    const std::vector<Input::LineDefinition>& possible_lines)
{
  auto [parsed_lines, unparsed_lines] =
      read_matching_lines_in_section(reader, section, possible_lines);

  // In this function, encountering unparsed lines is an error, so construct a nice message.
  if (unparsed_lines.size() > 0)
  {
    std::stringstream out;
    out << "Read failed in section " << std::quoted(section) << '\n';
    for (const auto& unparsed : unparsed_lines)
    {
      out << "  " << std::quoted(unparsed) << '\n';
    }
    out << "Valid lines are:\n";
    std::for_each(possible_lines.begin(), possible_lines.end(),
        [&](const Input::LineDefinition& def)
        {
          def.print(out);
          out << '\n';
        });
    FOUR_C_THROW(out.str().c_str());
  }

  return parsed_lines;
}


std::pair<std::vector<Input::LineDefinition>, std::vector<std::string>>
Core::IO::DatFileUtils::read_matching_lines_in_section(Core::IO::DatFileReader& reader,
    const std::string& section, const std::vector<Input::LineDefinition>& possible_lines)
{
  std::vector<std::string> unparsed_lines;
  std::vector<Input::LineDefinition> parsed_lines;

  Input::LineDefinition::ReadContext context{.input_file = reader.my_inputfile_name()};

  const auto process_line = [&](const std::string& input_line)
  {
    for (const auto& definition : possible_lines)
    {
      std::stringstream l{input_line};

      // Make a copy that potentially gets filled by the Read.
      auto parsed_definition = definition;
      if (parsed_definition.read(l, context))
      {
        parsed_lines.emplace_back(std::move(parsed_definition));
        return;
      }
    }
    // None of the possible lines matched.
    unparsed_lines.emplace_back(input_line);
  };

  for (const auto& input_line : reader.lines_in_section("--" + section))
  {
    process_line(std::string(input_line));
  }

  return {parsed_lines, unparsed_lines};
}
FOUR_C_NAMESPACE_CLOSE

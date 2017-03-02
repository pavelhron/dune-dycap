// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARSER_HH
#define DUNE_DYCAP_PARSER_HH

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

// read parameters from file
struct Parser {
  Parser(Dune::ParameterTree & param_): param(param_)
  {}
  void operator()(std::string n) {
    Dune::ParameterTreeParser parser;
    parser.readINITree(n, param);
  }

  Dune::ParameterTree & param;
};

void addConfigFiles(Dune::ParameterTree & param)
{
  std::vector<std::string> inputs =  param.get<std::vector<std::string> >("inputs", std::vector<std::string>());
  for_each( inputs.begin(), inputs.end(),Parser(param));
}

#endif

// lmr: runs lmscript `_lmcode_` bytecode against the native statistics core (no R, no Python, no Docker).
//
//   lmr --lmcode-schema           compact summary of the API (context, functions, types)
//   lmr --lmcode-schema --json    the embedded schema bytes: what lmscript's `_lmcode_ cli="lmr"` reads and compiles against
//                              (or save them for `_lmcode_ schema="..."`); --schema is an older alias
//   lmr --lmcode-signatures       canonical registry signatures (diffed against `lmscript --lmcode-signatures`)
//   lmr --generate-schema      the schema JSON derived from the registry (maintenance: scripts/update-schema.sh)
//   lmr --lmcode A.dslb [--lmcode B.dslb ...]  load and verify every file, then run them in slot order with shared host state
//   lmr --version
//
// Exit codes: 0 ok, 2 usage, 3 load/run error, 4 application failure (out.fail or failed out.check calls).
// stdout carries exactly one JSON result line; out.print text goes to stderr.
#include <iostream>
#include <map>
#include <optional>
#include <set>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "host.hpp"
#include "schema_embedded.hpp"

#ifndef LMR_VERSION
#define LMR_VERSION "0.5.2"
#endif

namespace {

using namespace rcli::host;

struct Parsed {
  std::string kind, type, name, returns;
  Strs params;
};

Parsed parse_signature(const std::string& s) {
  Parsed p;
  if (s.rfind("context ", 0) == 0) {
    p.kind = "context";
    size_t colon = s.find(':');
    p.name = s.substr(8, colon - 8);
    p.type = s.substr(colon + 2);
    return p;
  }
  bool is_fn = s.rfind("fn ", 0) == 0;
  p.kind = is_fn ? "fn" : "method";
  size_t start = is_fn ? 3 : 7;
  size_t open = s.find('('), close = s.find(')');
  std::string qual = s.substr(start, open - start);
  if (!is_fn) {
    size_t dot = qual.find('.');
    p.type = qual.substr(0, dot);
    p.name = qual.substr(dot + 1);
  } else {
    p.name = qual;
  }
  std::string inside = s.substr(open + 1, close - open - 1);
  size_t pos = 0;
  while (!inside.empty() && pos <= inside.size()) {
    size_t comma = inside.find(", ", pos);
    p.params.push_back(inside.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos));
    if (comma == std::string::npos) break;
    pos = comma + 2;
  }
  p.returns = s.substr(s.find("-> ") + 3);
  return p;
}

std::string overloads_json(const std::vector<Parsed>& list, const std::string& indent) {
  std::string out = "[";
  for (size_t i = 0; i < list.size(); ++i) {
    out += (i ? ",\n" + indent + " " : "") + std::string("{\"params\": [");
    for (size_t j = 0; j < list[i].params.size(); ++j) out += (j ? ", " : "") + dslvm::json_escape(list[i].params[j]);
    out += "], \"returns\": " + dslvm::json_escape(list[i].returns) + "}";
  }
  return out + "]";
}

// Deterministic, diff-friendly schema text: types, methods and functions in sorted order, one name per line.
std::string generate_schema(const dslvm::Registry& reg, const OperatorTable& ops) {
  std::map<std::string, std::map<std::string, std::vector<Parsed>>> types;
  std::map<std::string, std::vector<Parsed>> functions;
  std::map<std::string, std::string> context;
  for (const auto& s : reg.signatures()) {
    Parsed p = parse_signature(s);
    // The lmcode registry now supplies these standard StrList methods itself.
    // Keep the other standard declarations in the schema for compiler parity.
    if (p.kind == "method" && p.type == "StrList" &&
        (p.name == "len" || p.name == "at" || p.name == "contains" || p.name == "join")) continue;
    if (p.kind == "context") context[p.name] = p.type;
    else if (p.kind == "fn") functions[p.name].push_back(p);
    else types[p.type][p.name].push_back(p);
  }
  std::string out = "{\n  \"dsl.schema\": 2,\n  \"name\": \"lmr\",\n  \"context\": {";
  size_t i = 0;
  for (const auto& c : context) out += (i++ ? ", " : "") + dslvm::json_escape(c.first) + ": " + dslvm::json_escape(c.second);
  out += "},\n  \"types\": {\n";
  i = 0;
  for (const auto& t : types) {
    out += (i++ ? ",\n" : "") + std::string("    ") + dslvm::json_escape(t.first) + ": {\n";
    auto op = ops.find(t.first);
    if (op != ops.end()) {
      out += "      \"operators\": {";
      size_t k = 0;
      for (const auto& o : op->second) out += (k++ ? ", " : "") + dslvm::json_escape(o.first) + ": " + dslvm::json_escape(o.second);
      out += "},\n";
    }
    out += "      \"methods\": {\n";
    size_t j = 0;
    for (const auto& m : t.second)
      out += (j++ ? ",\n" : "") + std::string("        ") + dslvm::json_escape(m.first) + ": " + overloads_json(m.second, "        ");
    out += "\n      }\n    }";
  }
  out += "\n  },\n  \"functions\": {\n";
  i = 0;
  for (const auto& f : functions)
    out += (i++ ? ",\n" : "") + std::string("    ") + dslvm::json_escape(f.first) + ": " + overloads_json(f.second, "    ");
  out += "\n  }\n}\n";
  return out;
}

void print_summary(const dslvm::Registry& reg) {
  std::map<std::string, std::vector<std::string>> members;
  for (const auto& s : reg.signatures()) {
    if (s.rfind("context ", 0) == 0 || s.rfind("fn ", 0) == 0) { std::cout << s << "\n"; continue; }
    std::string body = s.substr(7);
    members[body.substr(0, body.find('.'))].push_back(body.substr(body.find('.') + 1));
  }
  for (const auto& t : members) {
    std::cout << t.first << ":\n";
    for (const auto& m : t.second) std::cout << "  " << m << "\n";
  }
}

std::string result_json(const Out& out) {
  std::string j = "{\"values\": {";
  for (size_t i = 0; i < out.values.size(); ++i)
    j += (i ? ", " : "") + dslvm::json_escape(out.values[i].first) + ": " + out.values[i].second;
  return j + "}, \"checks\": " + std::to_string(out.checks) + ", \"failures\": " + std::to_string(out.failures) + "}";
}

}  // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  // Byte-exact output: the schema is compared and hashed, and the JSON line must not gain a CR.
  _setmode(_fileno(stdout), _O_BINARY);
  _setmode(_fileno(stderr), _O_BINARY);
#endif
  Host host;
  dslvm::Registry reg;
  OperatorTable ops;
  Binder binder(reg, host, ops);
  try {
    bind_all(binder);
  } catch (const std::exception& e) {
    std::cerr << "lmr: registry error: " << e.what() << "\n";
    return 3;
  }
  const std::string usage =
      "usage: lmr --lmcode-schema [--json] | --lmcode-signatures | --generate-schema | --version | --lmcode FILE [FILE ...]\n";
  std::vector<std::string> args(argv + 1, argv + argc);
  if (args.size() == 1 && args[0] == "--version") {
    std::cout << "lmr " << LMR_VERSION << "\n";
    return 0;
  }
  if (!args.empty() && (args[0] == "--lmcode-schema" || args[0] == "--schema")) {
    if (args.size() == 1) { print_summary(reg); return std::cout ? 0 : 4; }
    if (args.size() == 2 && args[1] == "--json") {
      std::cout.write(reinterpret_cast<const char*>(rcli::embedded_schema), static_cast<std::streamsize>(rcli::embedded_schema_size));
      return std::cout ? 0 : 4;
    }
    std::cerr << usage;
    return 2;
  }
  if (args.size() == 1 && args[0] == "--lmcode-signatures") {
    for (const auto& s : reg.signatures()) std::cout << s << "\n";
    return 0;
  }
  if (args.size() == 1 && args[0] == "--generate-schema") {
    std::cout << generate_schema(reg, ops);
    return 0;
  }
  // Both argv shapes are legal: --lmcode FILE [--lmcode FILE ...] (builder slot fragments) and
  // --lmcode FILE [FILE ...] (plain multi-file runs). Slot order comes from the invocation metadata.
  std::vector<std::string> files;
  bool seen_flag = false;
  for (size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--lmcode") {
      seen_flag = true;
      continue;
    }
    if (args[i].empty() || args[i].rfind("--", 0) == 0) {
      std::cerr << usage;
      return dslvm::failure_frame("cli", dslvm::Error("CLI_USAGE", "invalid program file argument"));
    }
    if (!seen_flag) {
      std::cerr << usage;
      return dslvm::failure_frame("cli", dslvm::Error("CLI_USAGE", "expected --lmcode FILE [FILE ...]"));
    }
    files.push_back(args[i]);
  }
  if (files.empty()) {
    std::cerr << usage;
    return dslvm::failure_frame("cli", dslvm::Error("CLI_USAGE", "expected at least one program file"));
  }

  const std::string schema(reinterpret_cast<const char*>(rcli::embedded_schema), rcli::embedded_schema_size);
  // One invocation envelope (lmcode-cpp InvocationResults). Exit codes stay 0 ok, 2 cli (usage),
  // 3 load/run error, 4 application failure (failed out.check calls).
  std::optional<dslvm::InvocationResults> holder;
  try {
    holder.emplace(dslvm::InvocationResults::from_command(files));
  } catch (const dslvm::Error& e) {
    std::cerr << usage;
    return dslvm::failure_frame("cli", e);
  }
  auto& results = *holder;
  dslvm::ResultContract contract;
  try { contract = dslvm::ResultContract::from_schema(schema); }
  catch (const dslvm::Error& e) { results.fail(std::nullopt, "load", e); return results.write_frame(std::cout); }
  // Load and verify every chunk before anything runs: a broken later chunk must not
  // leave any earlier chunk's output behind.
  struct Part { size_t index; std::string file; dslvm::Program program; };
  std::vector<Part> programs;
  for (const auto& part : results.registrations()) {
    try {
      programs.push_back(Part{part.index, part.file, dslvm::Program::load(dslvm::read_file(part.file), schema, reg)});
    } catch (const dslvm::Error& e) {
      dslvm::Error located("LOAD_BAD_BYTECODE", part.file + ": " + e.what());
      located.code = e.code;
      std::cerr << located.message << "\n";
      results.fail(part.index, "load", located);
      return results.write_frame(std::cout);
    }
  }
  dslvm::Context context;
  context.set("out", &host.out);
  size_t failures_before = 0;
  for (auto& part : programs) {
    const size_t index = part.index;
    try {
      auto& sink = results.select(index, reg, contract);
      // A dry invocation validates and loads every chunk but runs none (stats runs write files).
      if (!results.dry()) part.program.run(context, reg);
      // Publish the run summary (values in first-set order, checks and failures) into the frame
      // result, with the same shape as the old one-line protocol's payload.
      std::string summary = result_json(host.out);
      dslvm::stdlib::json_detail::Parser parsed(summary);
      dslvm::stdlib::JsonNode node = parsed.val();
      parsed.ws();
      if (parsed.p != summary.size()) throw dslvm::Error("RUN_HOST", "result summary did not parse");
      sink.set(dslvm::Value::object("Json", &node));
      if (host.out.failures > failures_before) {
        // Application failure (failed out.check calls): the message text is unchanged from the
        // pre-envelope one-line protocol; counts are host-wide.
        dslvm::Error app("APP_FAIL", std::to_string(host.out.failures) + " of " + std::to_string(host.out.checks) +
                                     " checks failed; first: " + host.out.first_failure);
        std::cerr << app.message << "\n";
        results.fail(index, "app", app);
        return results.write_frame(std::cout);
      }
      results.complete(index);
    } catch (const dslvm::Error& e) {
      const char* phase = e.code == "APP_FAIL" ? "app" : "run";
      if (programs.size() > 1) {
        dslvm::Error located = e;
        located.message = part.file + ": " + e.message;
        std::cerr << located.message << "\n";
        results.fail(index, phase, located);
      } else {
        std::cerr << e.what() << "\n";
        results.fail(index, phase, e);
      }
      return results.write_frame(std::cout);
    }
    failures_before = host.out.failures;
  }
  return results.write_frame(std::cout);
}

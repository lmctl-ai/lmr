#include "host.hpp"

#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace rcli::host {

// ---- Frame and CSV (RFC 4180: quoted fields, doubled quotes, CRLF or LF) ----

const Column& Frame::col(const std::string& name) const {
  int i = find(name);
  if (i < 0) {
    std::string have;
    for (const auto& n : names) have += (have.empty() ? "" : ", ") + n;
    throw dslvm::Error("frame has no column '" + name + "' (columns: " + have + ")");
  }
  return cols[static_cast<size_t>(i)];
}

void Frame::set(const std::string& name, Column c) {
  if (!cols.empty() && c.size() != nrow())
    throw dslvm::Error("column '" + name + "' has " + std::to_string(c.size()) + " rows; the frame has " + std::to_string(nrow()));
  int i = find(name);
  if (i >= 0) cols[static_cast<size_t>(i)] = std::move(c);
  else { names.push_back(name); cols.push_back(std::move(c)); }
}

bool parse_number(const char* b, const char* e, double& out) {
  const char* p = b;
  if (p < e && *p == '-') ++p;
  if (p == e || !(std::isdigit(static_cast<unsigned char>(*p)) || *p == '.' || *p == 'i' || *p == 'I' || *p == 'n' || *p == 'N'))
    return false;
  if (e - p > 1 && p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) return false;
  const std::string s(b, e);  // strtod needs a terminator; the C locale is never changed
  char* end = nullptr;
  errno = 0;
  const double v = std::strtod(s.c_str(), &end);
  if (end != s.c_str() + s.size()) return false;
  if (errno == ERANGE && (std::isinf(v) || v == 0.0)) return false;
  out = v;
  return true;
}

static bool parse_double(const std::string& s, double& out) {
  if (s.empty()) return false;
  const char* b = s.data();
  if (*b == '+') ++b;
  return parse_number(b, s.data() + s.size(), out);
}

static bool is_na(const std::string& s) { return s.empty() || s == "NA" || s == "NaN"; }

static std::vector<Strs> parse_csv_records(const std::string& text, const std::string& path) {
  std::vector<Strs> rows;
  Strs row;
  std::string field;
  bool quoted = false, field_started = false, after_quote = false;
  size_t line = 1;
  auto end_field = [&] { row.push_back(field); field.clear(); field_started = false; after_quote = false; };
  auto end_row = [&] {
    end_field();
    if (!(row.size() == 1 && row[0].empty())) rows.push_back(row);  // skip blank lines
    row.clear();
  };
  for (size_t i = 0; i < text.size(); ++i) {
    char c = text[i];
    if (quoted) {
      if (c == '"') {
        if (i + 1 < text.size() && text[i + 1] == '"') { field += '"'; ++i; }
        else { quoted = false; after_quote = true; }
      } else {
        if (c == '\n') ++line;
        field += c;
      }
      continue;
    }
    if (c == '"') {
      if (field_started || after_quote) throw dslvm::Error("read_csv: " + path + " line " + std::to_string(line) + ": unexpected quote");
      quoted = true; field_started = true;
    } else if (c == ',') {
      end_field();
    } else if (c == '\n' || c == '\r') {
      if (c == '\r' && i + 1 < text.size() && text[i + 1] == '\n') ++i;
      end_row();
      ++line;
    } else {
      if (after_quote) throw dslvm::Error("read_csv: " + path + " line " + std::to_string(line) + ": text after a closing quote");
      field += c; field_started = true;
    }
  }
  if (quoted) throw dslvm::Error("read_csv: " + path + ": unterminated quoted field");
  if (field_started || !row.empty() || !field.empty()) end_row();
  return rows;
}

Frame read_csv(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw dslvm::Error("read_csv: cannot open " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  std::string text = ss.str();
  if (text.size() >= 3 && static_cast<unsigned char>(text[0]) == 0xEF && static_cast<unsigned char>(text[1]) == 0xBB &&
      static_cast<unsigned char>(text[2]) == 0xBF)
    text.erase(0, 3);
  auto rows = parse_csv_records(text, path);
  if (rows.empty()) throw dslvm::Error("read_csv: " + path + " has no header row");
  Frame f;
  const Strs& header = rows[0];
  for (size_t j = 0; j < header.size(); ++j) {
    if (header[j].empty()) throw dslvm::Error("read_csv: " + path + ": column " + std::to_string(j + 1) + " has an empty name");
    for (size_t k = 0; k < j; ++k)
      if (header[k] == header[j]) throw dslvm::Error("read_csv: " + path + ": duplicate column '" + header[j] + "'");
  }
  size_t n = rows.size() - 1;
  for (size_t r = 1; r < rows.size(); ++r)
    if (rows[r].size() != header.size())
      throw dslvm::Error("read_csv: " + path + ": record " + std::to_string(r + 1) + " has " + std::to_string(rows[r].size()) +
                         " fields; the header has " + std::to_string(header.size()));
  for (size_t j = 0; j < header.size(); ++j) {
    Column c;
    c.num.resize(n);
    bool numeric = true;
    for (size_t r = 0; r < n && numeric; ++r) {
      const std::string& s = rows[r + 1][j];
      double v;
      if (is_na(s)) c.num[r] = kNaN;
      else if (parse_double(s, v)) c.num[r] = v;
      else numeric = false;
    }
    if (!numeric) {
      c.text = true;
      c.num.clear();
      c.str.resize(n);
      for (size_t r = 0; r < n; ++r) c.str[r] = rows[r + 1][j];
    }
    f.names.push_back(header[j]);
    f.cols.push_back(std::move(c));
  }
  return f;
}

void check_output_path(const std::string& path, const char* what) {
  if (path.empty()) throw dslvm::Error(std::string(what) + ": output path must not be empty");
  if (path[0] == '/' || path[0] == '\\' || (path.size() > 1 && path[1] == ':'))
    throw dslvm::Error(std::string(what) + ": output path must be relative to the working directory: " + path);
  std::string part;
  std::stringstream ss(path);
  while (std::getline(ss, part, '/'))
    if (part == "..") throw dslvm::Error(std::string(what) + ": output path must not contain '..': " + path);
  if (path.back() == '/') throw dslvm::Error(std::string(what) + ": output path names a directory: " + path);
}

static std::string csv_quote(const std::string& s) {
  if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) { if (c == '"') out += '"'; out += c; }
  return out + "\"";
}

void write_csv(const Frame& f, const std::string& path) {
  check_output_path(path, "write_csv");
  std::string tmp = path + ".tmp";
  {
    std::ofstream o(tmp, std::ios::binary | std::ios::trunc);
    if (!o) throw dslvm::Error("write_csv: cannot write " + path);
    for (size_t j = 0; j < f.names.size(); ++j) o << (j ? "," : "") << csv_quote(f.names[j]);
    o << "\n";
    char buf[40];
    for (size_t r = 0; r < f.nrow(); ++r) {
      for (size_t j = 0; j < f.cols.size(); ++j) {
        if (j) o << ",";
        const Column& c = f.cols[j];
        if (c.text) o << csv_quote(c.str[r]);
        else if (std::isnan(c.num[r])) o << "NA";
        else { std::snprintf(buf, sizeof buf, "%.17g", c.num[r]); o << buf; }
      }
      o << "\n";
    }
    if (!o) throw dslvm::Error("write_csv: write failed for " + path);
  }
  if (std::rename(tmp.c_str(), path.c_str()) != 0) {
    std::remove(tmp.c_str());
    throw dslvm::Error("write_csv: cannot replace " + path);
  }
}

// ---- Driver ----

bool Driver::resume() {
  std::unique_lock<std::mutex> lock(mutex_);
  if (finished_) return false;
  if (!started_) {
    if (!algorithm_) throw dslvm::Error("internal: driver has no algorithm");
    started_ = true;
    worker_turn_ = true;
    worker_ = std::thread([this] {
      try {
        algorithm_();
      } catch (const Cancelled&) {
      } catch (...) {
        error_ = std::current_exception();
      }
      std::lock_guard<std::mutex> g(mutex_);
      finished_ = true;
      worker_turn_ = false;
      cv_.notify_all();
    });
  } else {
    worker_turn_ = true;
    cv_.notify_all();
  }
  cv_.wait(lock, [this] { return !worker_turn_; });
  if (finished_ && error_) {
    auto e = error_;
    error_ = nullptr;
    std::rethrow_exception(e);
  }
  return !finished_;
}

void Driver::yield() {
  std::unique_lock<std::mutex> lock(mutex_);
  if (cancel_) throw Cancelled{};
  worker_turn_ = false;
  cv_.notify_all();
  cv_.wait(lock, [this] { return worker_turn_; });
  if (cancel_) throw Cancelled{};
}

Driver::~Driver() {
  if (!worker_.joinable()) return;
  {
    std::lock_guard<std::mutex> g(mutex_);
    if (!finished_) { cancel_ = true; worker_turn_ = true; }
    cv_.notify_all();
  }
  worker_.join();
}

// ---- Out ----

void Out::set(const std::string& key, std::string json) {
  for (auto& kv : values)
    if (kv.first == key) { kv.second = std::move(json); return; }
  values.emplace_back(key, std::move(json));
}

// ---- Binder ----

static bool kind_ok(const Value& v, const std::string& t) {
  if (t == "float") return v.kind == Value::FLOAT;
  if (t == "int") return v.kind == Value::INT;
  if (t == "string") return v.kind == Value::STR;
  if (t == "bool") return v.kind == Value::BOOL;
  return v.kind == Value::OBJ && v.type == t && v.obj != nullptr;
}

template <class Call>
static Value guarded(const std::string& what, const Strs& params, Args& a, Call&& call) {
  if (a.size() != params.size()) throw dslvm::Error(what + ": expected " + std::to_string(params.size()) + " arguments");
  for (size_t i = 0; i < params.size(); ++i)
    if (!kind_ok(a[i], params[i])) throw dslvm::Error(what + " argument " + std::to_string(i + 1) + " must be " + params[i]);
  try {
    return call();
  } catch (const dslvm::Error&) {
    throw;
  } catch (const std::exception& e) {
    throw dslvm::Error(e.what());
  }
}

void Binder::fn(const std::string& name, Strs params, const std::string& returns, dslvm::Function f) {
  reg.function(name, params, returns, [name, params, f](Args& a) { return guarded(name, params, a, [&] { return f(a); }); });
}

void Binder::m(const std::string& type, const std::string& name, Strs params, const std::string& returns, dslvm::Method f) {
  std::string what = type + "." + name;
  reg.method(type, name, params, returns, [what, params, f](void* s, Args& a) {
    if (!s) throw dslvm::Error(what + ": null receiver");
    return guarded(what, params, a, [&] { return f(s, a); });
  });
}

// ---- JSON and formatting ----

std::string json_number(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  // Shortest significant-digit count that round-trips, printed without an exponent when |v| is in [1e-5, 1e17):
  // %g switches to exponent form once the decimal exponent reaches the precision, so the precision is raised to
  // cover the integer digits (10 -> "10", not "1e+01"; 2598960 -> "2598960").
  for (int p = 1; p <= 17; ++p) {
    char shorter[40];
    std::snprintf(shorter, sizeof shorter, "%.*g", p, v);
    if (std::strtod(shorter, nullptr) != v) continue;
    double mag = std::fabs(v);
    if (mag >= 1e-5 && mag < 1e17) {
      int exp10 = static_cast<int>(std::floor(std::log10(mag)));
      if (exp10 + 1 > p) std::snprintf(shorter, sizeof shorter, "%.*g", exp10 + 1, v);
    }
    return shorter;
  }
  return buf;
}

std::string json_vec(const Vec& v) {
  std::string out = "[";
  for (size_t i = 0; i < v.size(); ++i) out += (i ? ", " : "") + json_number(v[i]);
  return out + "]";
}

// Accepts literal text with exactly one conversion %[flags][width][.precision](f|F|e|E|g|G|d|i|s) ('%%' is literal).
// The caller passes the value already converted for the conversion kind it validated.
static size_t find_conversion(const std::string& spec, char& conv) {
  size_t pos = std::string::npos;
  for (size_t i = 0; i < spec.size(); ++i) {
    if (spec[i] != '%') continue;
    if (i + 1 < spec.size() && spec[i + 1] == '%') { ++i; continue; }
    if (pos != std::string::npos) throw dslvm::Error("fmt: format must contain exactly one conversion: " + spec);
    size_t j = i + 1;
    while (j < spec.size() && std::string("-+ 0#").find(spec[j]) != std::string::npos) ++j;
    size_t digits = 0;
    while (j < spec.size() && std::isdigit(static_cast<unsigned char>(spec[j]))) { ++j; ++digits; }
    if (digits > 3) throw dslvm::Error("fmt: width too large: " + spec);
    if (j < spec.size() && spec[j] == '.') {
      ++j; digits = 0;
      while (j < spec.size() && std::isdigit(static_cast<unsigned char>(spec[j]))) { ++j; ++digits; }
      if (digits > 2) throw dslvm::Error("fmt: precision too large: " + spec);
    }
    if (j >= spec.size() || std::string("fFeEgGdis").find(spec[j]) == std::string::npos)
      throw dslvm::Error("fmt: unsupported conversion in " + spec + " (use %f %e %g %d %s with flags, width, precision)");
    conv = spec[j];
    pos = i;
    i = j;
  }
  if (pos == std::string::npos) throw dslvm::Error("fmt: format has no conversion: " + spec);
  return pos;
}

std::string format_number(const std::string& spec, double v) {
  char conv = 0;
  find_conversion(spec, conv);
  if (conv == 's') throw dslvm::Error("fmt: %s needs a string value: " + spec);
  std::vector<char> buf(512);
  int n;
  if (conv == 'd' || conv == 'i') {
    if (!std::isfinite(v) || std::fabs(v) > 9.2e18) throw dslvm::Error("fmt: %d needs a finite integral value");
    n = std::snprintf(buf.data(), buf.size(), spec.c_str(), static_cast<long long>(v));
  } else {
    n = std::snprintf(buf.data(), buf.size(), spec.c_str(), v);
  }
  if (n < 0 || static_cast<size_t>(n) >= buf.size()) throw dslvm::Error("fmt: result too long");
  return std::string(buf.data(), static_cast<size_t>(n));
}

std::string format_string(const std::string& spec, const std::string& s) {
  char conv = 0;
  find_conversion(spec, conv);
  if (conv != 's') throw dslvm::Error("fmt: a string value needs %s: " + spec);
  std::vector<char> buf(s.size() + spec.size() + 1100);
  int n = std::snprintf(buf.data(), buf.size(), spec.c_str(), s.c_str());
  if (n < 0 || static_cast<size_t>(n) >= buf.size()) throw dslvm::Error("fmt: result too long");
  return std::string(buf.data(), static_cast<size_t>(n));
}

void bind_all(Binder& b) {
  bind_core(b);
  bind_dist(b);
  bind_stats(b);
  bind_models(b);
  bind_smooth(b);
  bind_ts(b);
  bind_cluster(b);
  bind_optim(b);
}

}  // namespace rcli::host

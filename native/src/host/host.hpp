// The rcli DSL host: shared state, object arenas, type-checked registration and the reverse-communication driver.
#pragma once
#include <condition_variable>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "dslvm.hpp"
#include "rcli/common.hpp"

namespace rcli::host {

using dslvm::Args;
using dslvm::Value;
using Strs = std::vector<std::string>;

// A table of named columns, each numeric (NaN = NA) or text.
struct Column {
  bool text = false;
  Vec num;
  Strs str;
  size_t size() const { return text ? str.size() : num.size(); }
};
struct Frame {
  Strs names;
  std::vector<Column> cols;
  size_t nrow() const { return cols.empty() ? 0 : cols[0].size(); }
  int find(const std::string& name) const {
    for (size_t i = 0; i < names.size(); ++i) if (names[i] == name) return static_cast<int>(i);
    return -1;
  }
  const Column& col(const std::string& name) const;
  void set(const std::string& name, Column c);
};
Frame read_csv(const std::string& path);
void write_csv(const Frame& f, const std::string& path);
// Output paths must be relative, stay below the working directory and name a file (no "..", no absolute path).
void check_output_path(const std::string& path, const char* what);

// Runs a callback-style algorithm on a worker thread and lets the DSL script drive it one evaluation at a time.
// Only one side runs at any moment, so execution is deterministic. Destroying a Driver whose algorithm is still
// waiting for a value cancels it: the pending callback throws Cancelled on the worker, which unwinds and is joined.
class Driver {
 public:
  struct Cancelled {};
  Driver() = default;
  Driver(const Driver&) = delete;
  Driver& operator=(const Driver&) = delete;
  ~Driver();
  // Starts `algorithm` on first resume(). The algorithm calls yield() at each evaluation request.
  void set_algorithm(std::function<void()> algorithm) { algorithm_ = std::move(algorithm); }
  bool started() const { return started_; }
  bool finished() const { return finished_; }
  // Script side: run the worker until it requests a value (returns true) or finishes (returns false; rethrows its error).
  bool resume();
  // Worker side: hand control to the script and wait to be resumed. Throws Cancelled when the driver is destroyed.
  void yield();

 private:
  std::function<void()> algorithm_;
  std::thread worker_;
  std::mutex mutex_;
  std::condition_variable cv_;
  bool started_ = false, finished_ = false, worker_turn_ = false, cancel_ = false;
  std::exception_ptr error_;
};

// Values published by out.set(...) into the final JSON result line.
struct Out {
  std::vector<std::pair<std::string, std::string>> values;  // key -> JSON text, in first-set order
  int checks = 0, failures = 0;
  std::string first_failure;
  void set(const std::string& key, std::string json);
};

struct Host {
  Out out;
  std::vector<std::shared_ptr<void>> objects;  // every object handed to a script lives until the host exits
  std::map<std::string, std::shared_ptr<void>> kept;  // keep()/kept_*() across chunks: name -> object
  std::map<std::string, std::string> kept_type;
  size_t object_limit = 2000000;

  template <class T>
  Value hold(const char* type, T value) {
    if (objects.size() >= object_limit) throw dslvm::Error(std::string("object limit exceeded (") + std::to_string(object_limit) + ")");
    auto p = std::make_shared<T>(std::move(value));
    objects.push_back(p);
    return Value::object(type, p.get());
  }
  // For types that are not movable (drivers): construct in place.
  template <class T, class... A>
  std::pair<Value, T*> emplace(const char* type, A&&... args) {
    if (objects.size() >= object_limit) throw dslvm::Error("object limit exceeded");
    auto p = std::make_shared<T>(std::forward<A>(args)...);
    objects.push_back(p);
    return {Value::object(type, p.get()), p.get()};
  }
  std::shared_ptr<void> owner(void* p) const {
    for (auto it = objects.rbegin(); it != objects.rend(); ++it) if (it->get() == p) return *it;
    return nullptr;
  }
};

// Operator mappings that the schema declares for a type (schema 2 "operators").
using OperatorTable = std::map<std::string, std::map<std::string, std::string>>;

// Type-checked registration: every argument is checked against the declared type before the callback runs
// (the VM trusts bytecode; this keeps a malformed or mismatched chunk from reaching host code), and rcli::Error or
// std::exception from numeric code become dslvm::Error (RUN_HOST) with the message unchanged.
class Binder {
 public:
  Binder(dslvm::Registry& r, Host& h, OperatorTable& ops) : reg(r), host(h), operators(ops) {}
  dslvm::Registry& reg;
  Host& host;
  OperatorTable& operators;
  void fn(const std::string& name, Strs params, const std::string& returns, dslvm::Function f);
  void m(const std::string& type, const std::string& name, Strs params, const std::string& returns, dslvm::Method f);
  void op(const std::string& type, const std::string& symbol, const std::string& method) { operators[type][symbol] = method; }
};

template <class T>
T& self(void* p) { return *static_cast<T*>(p); }
template <class T>
T& obj(const Value& v) { return *static_cast<T*>(v.obj); }

inline int64_t as_index(int64_t i, size_t n, const char* what) {
  if (i < 0 || static_cast<uint64_t>(i) >= n)
    throw dslvm::Error(std::string(what) + ": index " + std::to_string(i) + " out of range 0.." + std::to_string(static_cast<int64_t>(n) - 1));
  return i;
}
std::string json_number(double v);
std::string json_vec(const Vec& v);
std::string format_number(const std::string& spec, double v);         // one printf conversion, validated
std::string format_string(const std::string& spec, const std::string& s);  // one %s conversion

// Registration entry points, one per domain (bind_*.cpp).
void bind_core(Binder& b);     // Out, Vec, StrList, Frame, math, formatting, keep/kept
void bind_dist(Binder& b);     // Dist, Rng, special functions
void bind_stats(Binder& b);    // quantiles, hypothesis tests, FFT
void bind_models(Binder& b);   // lm, glm, nls
void bind_smooth(Binder& b);   // splines, approx, smooth.spline, loess
void bind_ts(Binder& b);       // acf, pacf, arima, arima_sim, stl
void bind_cluster(Binder& b);  // kmeans, hclust
void bind_optim(Binder& b);    // minimize, optimize, uniroot, integrate
void bind_all(Binder& b);

}  // namespace rcli::host

#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

namespace tomo {

class benchmark {
   public:
    using TTimePoint =
        std::chrono::time_point<std::chrono::high_resolution_clock>;

    benchmark(std::string title) : title_(title) {
        start_ = std::chrono::high_resolution_clock::now();
    }

    ~benchmark() {
        if (!silent_ && !finished_) {
            finish();
        }
    }

    void phase(std::string split_title) {
        split_title.resize(30, ' ');
        auto now = std::chrono::high_resolution_clock::now();
        splits_.push_back(make_pair(split_title, now));
    }

    void silence() { silent_ = true; }

    void finish() {
        if (silent_) return;

        finished_ = true;

        auto end = std::chrono::high_resolution_clock::now();
        auto total_ms =
            std::chrono::duration<double, std::milli>(end - start_).count();

        std::stringstream splitOutput;
        if (!splits_.empty()) {
            splits_.push_back(make_pair("", end));
            auto hline =
                "----------------------------------------------------------";
            splitOutput << "\n" << hline << '\n';
            for (unsigned int i = 0; i < splits_.size() - 1; ++i) {
                auto splitTime = splits_[i + 1].second - splits_[i].second;
                auto ms = std::chrono::duration<double, std::milli>(splitTime)
                              .count();
                splitOutput << std::fixed << std::setprecision(2)
                            << splits_[i].first << " \t" << ms << " ms"
                            << " \t" << (ms / total_ms) * 100 << "%\n";
            }
            splitOutput << hline;
        }

        std::cout << title_ << " total runtime: " << total_ms << " ms"
                       << splitOutput.str() << "\n";
    }

    auto splits() const { return splits_; }
    auto start() const { return start_; }

   private:
    std::vector<std::pair<std::string, TTimePoint>> splits_;
    std::string title_;
    bool silent_ = false;
    bool finished_ = false;
    TTimePoint start_;
};

} // namespace tomo

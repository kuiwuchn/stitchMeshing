#pragma once

#include "common.h"
#include <chrono>

template <typename TimeT = std::chrono::milliseconds> class Timer {
public:
    Timer() {
        start = std::chrono::system_clock::now();
    }

    size_t value() const {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        return (size_t) duration.count();
    }

    size_t reset() {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        start = now;
        return (size_t) duration.count();
    }

    void beginStage(const std::string &name) {
        reset();
        std::cout << name << " .. ";
        std::cout.flush();
    }

    void endStage(const std::string &str = "") {
        std::cout << "done. (took " << value() << " ms";
        if (!str.empty())
            std::cout << ", " << str;
        std::cout << ")" << std::endl;
    }
private:
    std::chrono::system_clock::time_point start;
};

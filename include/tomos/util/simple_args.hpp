#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace tomo {

struct options {
    int argc;
    char** argv;

    bool passed(std::string flag) {
        return std::find(argv, argv + argc, flag) !=
               (argv + argc);
    }

    std::string arg(std::string flag) {
        auto pos = std::find(argv, argv + argc, flag);
        if (pos == argv + argc || pos + 1 == argv + argc) {
            return "";
        }
        pos++;

        return std::string(*pos);
    }

    std::vector<std::string> args(std::string flag) {
        auto pos = std::find(argv, argv + argc, flag);
        std::vector<std::string> result;
        if (pos == argv + argc || pos + 1 == argv + argc) {
            return result;
        }
        pos++;
        while(pos != argv + argc && *pos[0] != '-') {
            result.push_back(std::string(*pos));
            pos++;
        }

        return result;
    }

    template <typename T>
    T arg_as(std::string flag) {
        auto value = std::stringstream(arg(flag));
        T x = {};
        value >> x;
        return x;
    }

    template <typename T>
    T arg_as_or(std::string flag, T alt) {
        if (!passed(flag)) {
            return alt;
        }
        auto value = std::stringstream(arg(flag));
        T x = {};
        value >> x;
        return x;
    }

    bool required_arguments(const std::vector<std::string>& args) {
        for (auto& arg : args) {
            if (!passed(arg)) {
                return false;
            }
        }
        return true;
    }
};

} // namespace tomo

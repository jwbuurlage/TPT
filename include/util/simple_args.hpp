#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace tomo {

struct options {
    int argc;
    char** argv;
};

bool passed(options opts, std::string flag) {
    return std::find(opts.argv, opts.argv + opts.argc, flag) !=
           (opts.argv + opts.argc);
}

std::string arg(options opts, std::string flag) {
    auto pos = std::find(opts.argv, opts.argv + opts.argc, flag);
    if (pos == opts.argv + opts.argc || pos + 1 == opts.argv + opts.argc) {
        return "";
    }
    pos++;

    return std::string(*pos);
}

template <typename T>
T arg_as(options opts, std::string flag) {
    auto value = std::stringstream(arg(opts, flag));
    T x = {};
    value >> x;
    return x;
}

bool required_arguments(options opts, std::vector<std::string> args) {
    for (auto& arg : args) {
        if (!passed(opts, arg)) {
            return false;
        }
    }
    return true;
}

} // namespace tomo

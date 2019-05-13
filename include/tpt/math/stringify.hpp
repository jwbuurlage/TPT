#include "../math/vector.hpp"
#include "../math.hpp"

#include <glm/gtx/string_cast.hpp>

namespace tpt {
namespace math {

template <dimension D, typename T>
auto to_string(vec<D, T> x) {
    return glm::to_string(x);
}

template <dimension D, typename T>
auto line_to_string(line<D, T> x) {
    return to_string(x.origin) + " + d. " + to_string(x.delta);
}

} // namespace math
} // namespace tpt

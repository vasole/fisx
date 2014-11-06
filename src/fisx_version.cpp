#include "fisx_version.h"

namespace fisx
{

const std::string & fisxVersion(void)
{
    static std::string result = FISX_VERSION_STR;
    return result;
}

} // namespace fisx

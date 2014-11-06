#ifndef FISX_VERSION_H
#define FISX_VERSION_H
#include <string>

#define FISX_VERSION_STR "1.0.0"
/*
   FISX_VERSION is (major << 16) + (minor << 8) + patch.
*/
#define FISX_VERSION_HEX 0x010000

namespace fisx
{

const std::string & fisxVersion(void);

} // namespace fisx


#endif //FISX_VERSION_H

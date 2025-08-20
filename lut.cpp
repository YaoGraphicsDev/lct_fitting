#include "lut.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

int LUT::write_to_exr(const std::string& name) {
    const char* err = nullptr;

    int ret = SaveEXR(reinterpret_cast<const float*>(tab.data()),
        static_cast<int>(width), static_cast<int>(height), 4/*rgba*/, 0, name.c_str(), &err);

    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            fprintf(stderr, "SaveEXR error: %s\n", err);
            FreeEXRErrorMessage(err);
        }
        return -1;
    }

    return 0;
}


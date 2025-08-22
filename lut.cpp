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

int LUT::read_from_exr(const std::string& name) {
    const char* err = nullptr;

    float* rgba = nullptr;
    int width, height;

    // Load the EXR file
    int ret = LoadEXR(&rgba, &width, &height, name.c_str(), &err);

    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            std::cout << "LoadEXR error: " << err << std::endl;
            FreeEXRErrorMessage(err);
        }
        return -1;
    }

    this->width = static_cast<size_t>(width);
    this->height = static_cast<size_t>(height);
    tab.resize(this->width * this->height * 4); // Assuming 4 channels (RGBA)

    std::memcpy(tab.data(), rgba, tab.size() * sizeof(float));

    // Free the memory allocated by TinyEXR
    free(rgba);

    return 0;
}


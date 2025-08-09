#pragma once

#include <vector>
#include <string>

#include "glm/glm.hpp"

struct LUT {
    LUT(size_t w, size_t h) : width(w), height(h), image(w* h, glm::vec4(0.0f)) {}

    glm::vec4& at(size_t row, size_t col) {
        assert(row < height && col < width);
        return image[row * width + col];
    }

	int write_to_exr(const std::string& name);

	size_t width;
	size_t height;
	std::vector<glm::vec4> image;
};

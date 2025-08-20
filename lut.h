#pragma once

#include <vector>
#include <string>
#include <iostream>

#include "glm/glm.hpp"

struct LUT {
    LUT(size_t w, size_t h) : width(w), height(h), tab(w* h, glm::vec4(0.0f)) {}

    glm::vec4& at(size_t row, size_t col) {
        assert(row < height && col < width);
        return tab[row * width + col];
    }

	int write_to_exr(const std::string& name);

	void print_element(int n, std::ostream& os) {
		for (size_t row = 0; row < height; ++row) {
			size_t offset = row * width;
			for (size_t col = 0; col < width; ++col) {
				os << tab[offset + col][n] << "\t";
			}
			os << std::endl;
		}
	}

	size_t width;
	size_t height;
	std::vector<glm::vec4> tab;
};

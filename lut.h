#pragma once

#include <vector>
#include <string>
#include <iostream>

#include "glm/glm.hpp"

struct LUT {
    LUT(size_t w, size_t h) : width(w), height(h), tab(w* h, glm::vec4(0.0f)) {}

	LUT(const std::string& name) {
		int err = read_from_exr(name);
		assert(err == 0);
		if (err != 0) {
			exit(1);
		}
	}

    glm::vec4& at(size_t row, size_t col) {
        assert(row < height && col < width);
        return tab[row * width + col];
    }

	int write_to_exr(const std::string& name);

	int read_from_exr(const std::string& name);

	void print_element(int n, std::ostream& os) {
		for (size_t row = 0; row < height; ++row) {
			size_t offset = row * width;
			for (size_t col = 0; col < width; ++col) {
				os << tab[offset + col][n] << "\t";
			}
			os << std::endl;
		}
	}

	glm::vec4 sample(glm::vec2 uv) {
		uv = glm::clamp(uv, glm::vec2(0.0f), glm::vec2(1.0f));
		glm::vec2 image_coord = uv * glm::vec2(width - 1, height - 1);
		glm::vec2 coord_floor = glm::floor(image_coord);
		int l = (int)coord_floor.x;
		int r = l + 1;
		int t = (int)coord_floor.y;
		int b = t + 1;

		glm::vec4 p00 = at(t, l);
		glm::vec4 p10 = at(t, r);
		glm::vec4 p01 = at(b, l);
		glm::vec4 p11 = at(b, r);
		float du = image_coord.x - l;
		float dv = image_coord.y - t;
		return glm::mix(glm::mix(p00, p10, du), glm::mix(p01, p11, du), dv);
	}

	size_t width;
	size_t height;
	std::vector<glm::vec4> tab;
};

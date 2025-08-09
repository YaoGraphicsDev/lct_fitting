#pragma once

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtc/random.hpp"
#include "glm/gtc/constants.hpp"
#include <random>

inline float pow2(float x) {
	return x * x;
}

// [0, 1]
inline float rand_float() {
	return (float)std::rand() / RAND_MAX;
}

inline glm::vec2 rand_vec2() {
	return glm::vec2(rand_float(), rand_float());
}

struct Stratified2DRNG {
	Stratified2DRNG(int strata) :
		strata(strata),
		strata2(strata * strata),
		cur_strata(0),
		strata_dim(1.0f / strata) {}

	glm::vec2 rand() {
		int row = cur_strata / strata;
		int col = cur_strata % strata;
		cur_strata++;
		cur_strata %= strata2;
		
		return (glm::vec2(col, row) + rand_vec2()) * strata_dim;
	}

	int strata;
	int strata2;
	int cur_strata;
	float strata_dim;
};



inline glm::vec3 random_on_hemisphere(glm::vec3 n) {
	glm::vec3 v = glm::sphericalRand(1.0f);
	if (glm::dot(n, v) < 0.0f) {
		return -v;
	}
	else {
		return v;
	}
}

inline glm::vec2 sample_unit_disk(glm::vec2 u) {
	float r = std::sqrt(u[0]);
	float theta = 2 * glm::pi<float>() * u[1];
	return glm::vec2(r * glm::cos(theta), r * std::sin(theta));
}

// malley's method
inline glm::vec3 sample_cosine_hemisphere(glm::vec2 u) {
	glm::vec2 d = sample_unit_disk(u);
	float z = glm::sqrt(glm::max(0.0f, 1 - d.x * d.x - d.y * d.y));
	return glm::vec3(d.x, d.y, z);
}

inline glm::mat3 build_tbn_z(glm::vec3 n) {
	glm::vec3 z(0.0f, 0.0f, 1.0f);
	glm::vec3 x(1.0f, 0.0f, 0.0f);
	glm::vec3 t = glm::cross(z, n);
	float t_len = glm::length(t);
	if (t_len < std::numeric_limits<float>::epsilon()) {
		t = x;
	}
	else {
		t = t / t_len;
	}
	glm::vec3 b = glm::cross(n, t);
	return glm::mat3(t, b, n);
}

// Unit Vector Arithmic
class UVA {
public:
	inline static float cos_theta(glm::vec3 v) {
		return v.z;
	}

	inline static float cos_theta_abs(glm::vec3 v) {
		return std::abs(cos_theta(v));
	}

	inline static float cos2_theta(glm::vec3 v) {
		return glm::clamp(cos_theta(v) * cos_theta(v), 0.0f, 1.0f);
	}

	inline static float sin_theta(glm::vec3 v) {
		return glm::sqrt(sin2_theta(v));
	}

	inline static float sin2_theta(glm::vec3 v) {
		return glm::clamp(1.0f - cos2_theta(v), 0.0f, 1.0f);
	}

	inline static float tan2_theta(glm::vec3 v) {
		return sin2_theta(v) / cos2_theta(v);
	}

	inline static float cos_phi(glm::vec3 v) {
		float s = sin_theta(v);
		return (s == 0.0f) ? 1.0f : glm::clamp(v.x / s, -1.0f, 1.0f);
	}

	inline static float sin_phi(glm::vec3 v) {
		float s = sin_theta(v);
		return (s == 0.0f) ? 0.0f : glm::clamp(v.y / s, -1.0f, 1.0f);
	}

	inline static float cos2_phi(glm::vec3 v) {
		float c = cos_phi(v);
		return c * c;
	}

	inline static float sin2_phi(glm::vec3 v) {
		float s = sin_phi(v);
		return s * s;
	}

	inline static float cos(glm::vec3 v1, glm::vec3 v2) {
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	inline static bool same_hemisphere(glm::vec3 v1, glm::vec3 v2) {
		return v1.z * v2.z > 0.0f;
	}

	inline static glm::vec3 spherical(float theta, float phi) {
		float sin_theta = glm::sin(theta);
		return glm::normalize(glm::vec3(sin_theta * glm::cos(phi), sin_theta * glm::sin(phi), glm::cos(theta)));
	}
};
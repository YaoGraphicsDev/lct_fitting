#pragma once

#include "utils.h"

#include "glm/glm.hpp"
#include "glm/gtc/constants.hpp"

#include <vector>

class LCT {
public:
	LCT(glm::vec3 X, glm::vec3 Y, glm::vec3 Z,
		float m11, float m13, float m22, float albedo) {

		// scaling and shearing. column major
		m[0] = { m11, 0.0f, 0.0f };
		m[1] = { 0.0f, m22, 0.0f };
		m[2] = { m13, 0.0f, 1.0f };
		
		m_params = { m11, m13, m22 };
		
		// rotation
		rotation = glm::mat3(X, Y, Z);
		m = rotation * m;

		inv_m = glm::inverse(m);
		det_inv_m = 1.0f / glm::abs(glm::determinant(m));
		_albedo = albedo;
	};

	LCT(glm::mat3 R, Vector<float, 3> ms, float albedo) 
		: LCT(R[0], R[1], R[2], ms[0], ms[1], ms[2], albedo) {}

	// This ctor does not properly set m_params and rotation
	LCT(Vector<float, 5> p, float albedo) {
		m[0] = { p[0], 0.0f, p[3]};
		m[1] = { 0.0f, p[2], 0.0f};
		m[2] = { p[1], 0.0f, p[4]};

		m_params = { 1.0f, 0.0f, 1.0f };
		rotation = glm::mat3(1.0f);

		inv_m = glm::inverse(m);
		det_inv_m = 1.0f / glm::abs(glm::determinant(m));
		_albedo = albedo;
	}

	struct EvalResult {
		float pdf;
		float brdfxcos_l;
	};
	EvalResult eval(glm::vec3 l) {
		glm::vec3 lo = inv_m * l;
		float inv_lo_len = 1.0f / glm::length(lo);
		lo *= inv_lo_len;
		float d_o = glm::one_over_pi<float>() * glm::clamp(UVA::cos_theta(lo), 0.0f, 1.0f);
		float J = det_inv_m * inv_lo_len * inv_lo_len * inv_lo_len;
		float pdf = d_o * J;
		assert(pdf >= 0.0f);
		return {pdf, pdf * _albedo};
	}

	glm::vec3 sample(glm::vec2 u) {
		glm::vec3 l = sample_cosine_hemisphere(u);
		//if ((m * l).z < 0.0f) {
		//	std::cout << "(ml).z = " << (m*l).z << std::endl;
		//}
		
		return glm::normalize(m * l);
	}

	// for debug only
	float albedo_sphere(int n = 64) {
		float sum = 0.0f;
		float d_theta = glm::pi<float>() / n;
		float d_phi = glm::two_pi<float>() / n;
		for (int i = 0; i < n; ++i) {
			float theta = (i + 0.5f) * d_theta;
			float cos_theta = glm::cos(theta);
			float sin_theta = glm::sin(theta);
			float d = sin_theta * d_phi * d_theta;
			for (int j = 0; j < n; ++j) {
				float phi = (j + 0.5f) * d_phi;
				glm::vec3 l(sin_theta * glm::cos(phi), sin_theta * glm::sin(phi), cos_theta);
				sum += eval(l).brdfxcos_l * d;
			}
		}

		return sum;
	}

	float albedo_hemisphere(int n = 64) {
		float sum = 0.0f;
		float d_theta = glm::half_pi<float>() / n;
		float d_phi = glm::two_pi<float>() / n;
		for (int i = 0; i < n; ++i) {
			float theta = (i + 0.5f) * d_theta;
			float cos_theta = glm::cos(theta);
			float sin_theta = glm::sin(theta);
			float d = sin_theta * d_phi * d_theta;
			for (int j = 0; j < n; ++j) {
				float phi = (j + 0.5f) * d_phi;
				glm::vec3 l(sin_theta * glm::cos(phi), sin_theta * glm::sin(phi), cos_theta);
				sum += eval(l).brdfxcos_l * d;
			}
		}

		return sum;
	}

// private:
	Vector<float, 3> m_params;
	glm::mat3 rotation;
	glm::mat3 m;
	glm::mat3 inv_m;
	float det_inv_m;
	float _albedo;
};
#pragma once

#include "utils.h"

#include <iostream>
#include <array>

class GGXDistribution {
public:
	GGXDistribution(float a) : alpha_x(glm::max(a, 0.0001f)), alpha_y(glm::max(a, 0.0001f)) {
		GGXDistribution(a, a);
	}

	GGXDistribution(float ax, float ay) : alpha_x(glm::max(ax, 0.0001f)), alpha_y(glm::max(ay, 0.0001f)) {}

	inline float D(glm::vec3 wm) {
		float k = glm::pi<float>() * alpha_x * alpha_y;
		float cos2_theta = UVA::cos2_theta(wm);
		float sin2_theta = UVA::sin2_theta(wm);
		float cos_phi = UVA::cos_phi(wm);
		float sin_phi = UVA::sin_phi(wm);
		float e = cos2_theta + sin2_theta * (pow2(cos_phi / alpha_x) + pow2(sin_phi / alpha_y));
		float result = 1.0f / (k * pow2(e));
		return result;
	}

	inline float D_visible(glm::vec3 w, glm::vec3 wm) {
		// safe way to compute G1/cos(theta)
		float cos = UVA::cos_theta_abs(w);
		float cos2 = UVA::cos2_theta(w);
		float sin2 = UVA::sin2_theta(w);
		float G1_over_cos =  2.0f / (cos + glm::sqrt(cos2 + alpha2(w) * sin2));
		float result = G1_over_cos * D(wm) * glm::max(0.0f, glm::dot(w, wm));
		return result;
	}

	inline float alpha2(glm::vec3 w) {
		return pow2(alpha_x) * UVA::cos2_phi(w) + pow2(alpha_y) * UVA::sin2_phi(w);
	}

	inline float lambda(glm::vec3 w) {
		float tan2 = UVA::tan2_theta(w);

		// perfectly grazing angle when cosine = 0
		// behaviour of a float devided by zero is implementation dependent
		// could be -NaN, NaN, inf, -inf.
		// return inf
		if (std::isnan(tan2)) {
			// expect lambda to be in denominator of G and G1, hence result in 0 value
			return std::numeric_limits<float>::infinity();
		}

		return (glm::sqrt(1 + alpha2(w) * tan2) - 1.0f) * 0.5f;
	}

	// conceptually, do not pass w that lies in xy plane
	inline float G1(glm::vec3 w) {
		return 1.0f / (1.0f + lambda(w));
	}
	
	inline float G(glm::vec3 wi, glm::vec3 wo) {
		return 1.0f / (1.0f + lambda(wi) + lambda(wo));
	}

	// v pointing outwards
	// sample normals
	inline glm::vec3 sample(glm::vec2 u, glm::vec3 v) {
		assert(alpha_x == alpha_y); // isotropic for the moment
		float alpha = alpha_x;
		float r = alpha * sqrtf(u.x / (1.0f - u.x));
		float phi = glm::two_pi<float>() * u.y;
		glm::vec3 m = glm::normalize(glm::vec3(r * glm::cos(phi), r * glm::sin(phi), 1.0f));
		glm::vec3 l = glm::reflect(-v, m);
		return glm::normalize(l);
	}

	//inline glm::vec3 sample_wm(glm::vec2 u, glm::vec3 w) {
	//	// "Sampling the GGX Distribution of Visible Normals" fig.4 & appendix A
	//	glm::vec3 Vh = glm::normalize(glm::vec3(w.x * alpha_x, w.y * alpha_y, w.z));
	//	float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
	//	float inv_lensq = glm::inversesqrt(lensq);
	//	
	//	// build projection space basis
	//	glm::vec3 T1 = lensq > 0.0f ? glm::vec3(-Vh.y, Vh.x, 0.0f) * inv_lensq : glm::vec3(1, 0, 0);
	//	glm::vec3 T2 = cross(Vh, T1);
	//	glm::mat3 basis = glm::mat3(T1, T2, Vh);
	//	glm::vec2 t1t2 = sample_unit_disk(u);
	//	float s = 0.5f * (1.0f + UVA::cos_theta(Vh));
	//	float t1 = t1t2.x;
	//	float t2 = t1t2.y;
	//	t2 = glm::mix(glm::sqrt(1.0f - t1 * t1), t2, s);
	//	
	//	glm::vec3 coord(t1, t2, glm::sqrt(glm::max(0.0f, 1.0f - t1 * t1 - t2 * t2)));
	//	glm::vec3 Nh = basis * coord;
	//	
	//	glm::vec3 Ne = glm::normalize(glm::vec3(Nh.x * alpha_x, Nh.y * alpha_y, std::max(0.0f, Nh.z)));
	//	return Ne;
	//}

public:
	float alpha_x;
	float alpha_y;
};


class BRDF {
public:
	BRDF(float alpha) : ggx(alpha) {};

	struct EvalResult {
		float pdf;
		float brdfxcos_l;
	};
	EvalResult eval(glm::vec3 v, glm::vec3 l) {
		glm::vec3 m = glm::normalize(v + l);
		if (v.z <= 0.0f || l.z < 0.0f || m.z <= 0.0f) {
			return { 0.0f, 0.0f };
		}

		float D = ggx.D(m);
		float pdf = D * UVA::cos_theta(m) / (4.0f * UVA::cos(v, m));
		assert(pdf >= 0.0f);
		
		float brdfxcos_l = D * ggx.G(v, l) / (4.0f * UVA::cos_theta(v));
		return { pdf, brdfxcos_l };
	}

	struct SampleResult {
		bool lv_same_hemisphere;
		glm::vec3 l;
	};
	SampleResult sample(glm::vec2 u, glm::vec3 v) {
		glm::vec3 l = ggx.sample(u, v);
		return { UVA::same_hemisphere(l, v), l };
	}

	float albedo_1(glm::vec3 v, int n = 64) {
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
				sum += eval(v, l).brdfxcos_l * d;
			}
		}

		return sum;
	}

	struct AlbedoContext {
		float albedo;
		glm::vec3 l_avg;
	};
	AlbedoContext albedo_2(glm::vec3 v, int n = 64) {
		float sum = 0.0f;
		glm::vec3 dir(0.0f);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				glm::vec2 u((i + 0.5f) / n, (j + 0.5f) / n);

				SampleResult s_result = sample(u, v);
				if (!s_result.lv_same_hemisphere) {
					sum += 0.0f;
					continue;
				}

				EvalResult e_result = eval(v, s_result.l);

				// accumulate
				assert(e_result.pdf > 0.0f);
				assert(e_result.brdfxcos_l >= 0.0f);

				sum += e_result.brdfxcos_l / e_result.pdf;
				dir += e_result.brdfxcos_l * s_result.l / e_result.pdf;
			}
		}

		// keep l in the same plane as v and +z
		glm::vec3 n_vz = glm::cross(v, glm::vec3(0.0f, 0.0f, 1.0f));
		if (glm::length(n_vz) < 0.0001) {
			dir = glm::vec3(0.0f, 0.0f, 1.0f);
		}
		else {
			n_vz = glm::normalize(n_vz);
			dir = dir - glm::dot(n_vz, dir) * n_vz;
		}
		
		return { sum / (float)(n * n), glm::normalize(dir) };
	}

	GGXDistribution ggx;
};

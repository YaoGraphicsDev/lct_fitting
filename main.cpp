// #include "nelder_mead.h"
#include "nelder_mead.hpp"
#include "brdf.h"
#include "glm/glm.hpp"
#include "lct.hpp"
#include "print.hpp"
#include "lut.h"

const int n_samples = 64;

float error_func(Vector<5> coord, glm::vec3 v, BRDF& brdf, float albedo) {
	// BRDF brdf(alpha);

	// float albedo = brdf.albedo_2(v);

	// std::cout << "coord = " << coord << std::endl;
    LCT lct(coord[0], coord[1], coord[2], coord[3], coord[4], albedo);

	float e3 = 0.0f;

	for (int i = 0; i < n_samples; ++i) {
		for (int j = 0; j < n_samples; ++j) {
			glm::vec2 u((i + 0.5f) / n_samples, (j + 0.5f) / n_samples);

			// importance sample LTC
			{
				// sample
				glm::vec3 l = lct.sample(u);

				// error with MIS weight
				BRDF::EvalResult eval_brdf = brdf.eval(v, l);
				LCT::EvalResult eval_lct = lct.eval(l);
				assert(eval_brdf.pdf >= 0.0f);
				assert(eval_lct.pdf >= 0.0f);
				assert(eval_brdf.pdf + eval_lct.pdf > 0.0f);
				assert(eval_brdf.brdfxcos_l >= 0.0f);
				assert(eval_lct.brdfxcos_l >= 0.0f);

				double e = fabsf(eval_brdf.brdfxcos_l - eval_lct.brdfxcos_l);
				e = e * e * e;
				e3 += e / (eval_lct.pdf + eval_brdf.pdf);
				//std::cout << "e1 = " << e / (eval_lct.pdf + eval_brdf.pdf) << std::endl;
			}

			// importance sample BRDF
			{
				// sample
				BRDF::SampleResult s_result = brdf.sample(u, v);
				if (!s_result.lv_same_hemisphere) {
					e3 += 0.0f;
					continue;
				}
				glm::vec3 l = s_result.l;

				// error with MIS weight
				BRDF::EvalResult eval_brdf = brdf.eval(v, l);
				LCT::EvalResult eval_lct = lct.eval(l);
				assert(eval_brdf.pdf >= 0.0f);
				assert(eval_lct.pdf >= 0.0f);
				assert(eval_brdf.pdf + eval_lct.pdf > 0.0f);
				assert(eval_brdf.brdfxcos_l >= 0.0f);
				assert(eval_lct.brdfxcos_l >= 0.0f);

				double e = fabsf(eval_brdf.brdfxcos_l - eval_lct.brdfxcos_l);
				e = e * e * e;
				e3 += e / (eval_lct.pdf + eval_brdf.pdf);
				// std::cout << "e2 = " << e / (eval_lct.pdf + eval_brdf.pdf) << std::endl;
			}
		}
	}

	e3 /= (float)(n_samples * n_samples);

	return e3;

}

void compare_lct_brdf(BRDF& brdf, LCT& lct, glm::vec3 v, int n = 32) {
	std::cout << "view direction v = " << v.x << ", " << v.y << ", " << v.z << std::endl;
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			float theta = glm::half_pi<float>() * glm::clamp((float)i / n, 0.0f, 1.0f);
			float phi = glm::two_pi<float>() * glm::clamp((float)j / n, 0.0f, 1.0f);
			glm::vec3 l = UVA::spherical(theta, phi);
			std::cout << "light direction l = " << l.x << ", " << l.y << ", " << l.z << std::endl;
			std::cout << "theta = " << glm::clamp((float)i / n, 0.0f, 1.0f) * 0.5f << " pi" << std::endl;
			std::cout << "phi = " << glm::clamp((float)j / n, 0.0f, 1.0f) * 2.0f << " pi" << std::endl;
			LCT::EvalResult lct_eval = lct.eval(l);
			BRDF::EvalResult brdf_eval = brdf.eval(v, l);
			std::cout << "lct:\t" << lct_eval.brdfxcos_l << std::endl;
			std::cout << "brdf:\t" << brdf_eval.brdfxcos_l << std::endl;
			std::cout << "diff = \t" << lct_eval.brdfxcos_l - brdf_eval.brdfxcos_l << std::endl;
		}
	}
}

struct FitResult {
	Vector<5> opti_x;
	float albedo;
};
FitResult fit(float theta_v, float alpha, Vector<5> start_x) {
	const float delta = 0.5f;
	const float tolerance = 0.00001;
	const int max_iter = 256;

	glm::vec3 v = UVA::spherical(theta_v, 0.0f);

	// align the center of cosine hemisphere starting point
	float sin_theta_v = glm::sin(theta_v);
	float cos_theta_v = v.z;
	
	// Vector<5> start_x = { cos_theta_v, sin_theta_v, 1.0f, -sin_theta_v, cos_theta_v };
	// Vector<5> start_x = { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f };
	BRDF brdf(alpha);
	float albedo = brdf.albedo_2(v);
	NelderMead<5>::Opti opti = NelderMead<5>::optimize(std::bind(error_func, std::placeholders::_1, v, brdf, albedo),
		start_x, delta, tolerance, max_iter);

	assert(opti.n_iter < max_iter);
	assert(opti.error < tolerance);

	return { opti.x, albedo };
}

LUT gen_lut(int resolution) {
	LUT lut(resolution, resolution);

	for (int i = 0; i < resolution; ++i) {
		float roughness = (float)i / (resolution - 1);
		float alpha = glm::clamp(roughness * roughness, 0.0001f, 1.0f);
		for (int j = 0; j < resolution; ++j) {
			float theta = (float)j / (resolution - 1) * glm::half_pi<float>();
			theta = glm::clamp(theta, 0.0f, glm::half_pi<float>() * 0.9999f); // theta = 0.5pi is a singularity
			// TODO: not done
		}
	}

	return lut;
}

int main() {
	float theta_v = glm::half_pi<float>() * 0.9f;
	glm::vec3 v = UVA::spherical(theta_v, 0.0f);
	float alpha = 0.001f;

	FitResult f_res = fit(theta_v, alpha, { 1.0f, 0.0f, 1.0f, 0.0f, 1.0f });
	LCT lct(f_res.opti_x, f_res.albedo);

	BRDF brdf(alpha);
	print_brdf(std::string(OUT_DATA_DIR) + "/brdf.dat", brdf, v);
	print_lct(std::string(OUT_DATA_DIR) + "/lct.dat", lct);

    return 0;
}
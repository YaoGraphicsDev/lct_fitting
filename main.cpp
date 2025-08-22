// #include "nelder_mead.h"
#include "nelder_mead.hpp"
#include "brdf.h"
#include "glm/glm.hpp"
#include "lct.hpp"
#include "print.hpp"
#include "lut.h"
#include "CLI.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/quaternion.hpp"

#include <tuple>

const int n_samples = 64;
const float delta = 0.05f;
const float tolerance = 1.0E-8;
const int max_iter = 128;

double error_func(glm::mat3 rotation, Vector<float, 3> coord, glm::vec3 v, BRDF& brdf, float albedo) {
	// BRDF brdf(alpha);

	// float albedo = brdf.albedo_2(v);

	// std::cout << "coord = " << coord << std::endl;
    LCT lct(rotation, coord, albedo);

	double e2 = 0.0f;

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
				e = e * e;
				e2 += e / (eval_lct.pdf + eval_brdf.pdf);
				//std::cout << "e1 = " << e / (eval_lct.pdf + eval_brdf.pdf) << std::endl;
			}

			// importance sample BRDF
			{
				// sample
				BRDF::SampleResult s_result = brdf.sample(u, v);
				if (!s_result.lv_same_hemisphere) {
					e2 += 0.0f;
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
				e = e * e;
				e2 += e / (eval_lct.pdf + eval_brdf.pdf);
				// std::cout << "e2 = " << e / (eval_lct.pdf + eval_brdf.pdf) << std::endl;
			}
		}
	}

	e2 /= (float)(n_samples * n_samples);

	return e2;

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
	glm::mat3 rotation;
	Vector<float, 3> opti_m_params;
	float albedo;
};
FitResult fit(glm::vec3 v, BRDF& brdf, float albedo, glm::vec3 l_avg, Vector<float, 3> start_x) {
	l_avg.y = 0.0f;
	std::cout << "ggx albedo = " << albedo << std::endl;
	std::cout << "ggx average direction = " << l_avg.x << ", " << l_avg.y << ", " << l_avg.z << std::endl;

	glm::vec3 t1(l_avg.z, 0.0f, -l_avg.x);
	glm::vec3 t2(0.0f, 1.0f, 0.0f);
	glm::mat3 rotation(t1, t2, l_avg);

	NelderMead<3>::Opti opti = NelderMead<3>::optimize(
		std::bind(error_func, rotation, std::placeholders::_1, v, brdf, albedo),
		start_x, delta, tolerance, max_iter);

	std::cout << "optimization result:" << std::endl;
	std::cout << "x = " << opti.x << std::endl;
	std::cout << "error = " << opti.error << std::endl;
	std::cout << "iterations = " << opti.n_iter << std::endl;

	assert(opti.n_iter < max_iter);
	assert(opti.error < tolerance);

	return { rotation, opti.x, albedo };
}

struct TabRange {
	int resolution;
	float min;
	float max;
};
std::vector<LUT> gen_lut(TabRange rough_range, TabRange theta_mult_range) {
	LUT lut_m(rough_range.resolution, theta_mult_range.resolution); // m11, m13, m22, m23
	
	LUT lut_b(rough_range.resolution, theta_mult_range.resolution);

	LUT lut_bm_0(rough_range.resolution, theta_mult_range.resolution);
	LUT lut_bm_1(rough_range.resolution, theta_mult_range.resolution);

	LUT lut_inv_0(rough_range.resolution, theta_mult_range.resolution);
	LUT lut_inv_1(rough_range.resolution, theta_mult_range.resolution);

	LUT lut_norm_inv(rough_range.resolution, theta_mult_range.resolution);

	Vector<float, 3> start_m_params;
	for (int r = rough_range.resolution - 1; r >= 0; --r) {
		float roughness = glm::clamp((float)r / (rough_range.resolution - 1), rough_range.min, rough_range.max);
		float alpha = roughness * roughness;
		BRDF brdf(alpha);
		for (int t = 0; t < theta_mult_range.resolution; ++t) {
			float theta_v_mult = glm::clamp((float)t / (theta_mult_range.resolution - 1), theta_mult_range.min, theta_mult_range.max);
			float theta_v = theta_v_mult * glm::half_pi<float>(); // theta = 0.5pi is a singularity
			glm::vec3 v = UVA::spherical(theta_v, 0.0f);
			BRDF::AlbedoContext a_ctx = brdf.albedo_2(v);

			std::cout << "roughness = " << roughness << ", theta = " << theta_v_mult << " * half pi" << std::endl;

			if (t == 0) {
				if (r == rough_range.resolution - 1) {
					start_m_params = { 1.0f, 0.0f, 1.0f };
				}
				else {
					start_m_params = {
						lut_m.at(t, r + 1)[0],
						lut_m.at(t, r + 1)[1],
						lut_m.at(t, r + 1)[2] };
				}
			}
			std::cout << "start_ms at " << start_m_params << std::endl;

			FitResult f_res = fit(v, brdf, a_ctx.albedo, a_ctx.l_avg, start_m_params);
			f_res.opti_m_params[0] = glm::abs(f_res.opti_m_params[0]); // prevent x and y from inverted scaling
			f_res.opti_m_params[2] = glm::abs(f_res.opti_m_params[2]);

			start_m_params = f_res.opti_m_params;

			// save the result of this round
			lut_m.at(t, r) = glm::vec4(
				f_res.opti_m_params[0],
				f_res.opti_m_params[1],
				f_res.opti_m_params[2], 0.0f);
			lut_b.at(t, r) = glm::vec4(
				f_res.rotation[0][2],
				f_res.rotation[2][2],
				0.0f, 0.0f);

			LCT lct(f_res.rotation, f_res.opti_m_params, a_ctx.albedo);
			glm::mat3 m = lct.m;
			lut_bm_0.at(t, r) = glm::vec4(m[0][0], m[2][0], m[1][1], 0.0f);
			lut_bm_1.at(t, r) = glm::vec4(m[0][2], m[2][2], a_ctx.albedo, 0.0f);
			
			glm::mat3 m_inv = glm::inverse(m);
			lut_inv_0.at(t, r) = glm::vec4(m_inv[0][0], m_inv[2][0], m_inv[1][1], 0.0f);
			lut_inv_1.at(t, r) = glm::vec4(m_inv[0][2], m_inv[2][2], a_ctx.albedo, 0.0f);

			glm::mat3 m_norm = m / m[2][2];
			float a = m_norm[0][0];
			float b = m_norm[0][2];
			float c = m_norm[1][1];
			float d = m_norm[2][0];
			lut_norm_inv.at(t, r) = glm::vec4(a, -b, (a - b * d) / c, -d);
		}
	}

	return { lut_m, lut_b, lut_bm_0, lut_bm_1, lut_inv_0, lut_inv_1, lut_norm_inv };
}

int main(int argc, char** argv) {
	CLI::App app;

	// output path
	std::string out_dir;
	app.add_option("-d,--output-dir", out_dir, "output directory")
		->check(CLI::ExistingDirectory)
		->default_val(OUT_DIR);

	auto dist_ggx_cmd = app.add_subcommand("dist-ggx", "Print values of GGX BRDF on upper hemisphere to a file");
	struct DistGGXParams {
		int resolution = 0;
		float theta_mult = 0.0f;
		float roughness = 0.0f;
		std::string out_filename = "";
	};
	DistGGXParams ggx_params;
	dist_ggx_cmd->add_option("-r,--res", ggx_params.resolution, "sampling resolution of theta and phi")
		->default_val(256);
	dist_ggx_cmd->add_option("-t,--theta-multiplier", ggx_params.theta_mult, "theta multiplier. Viewing angle theta_v = multiplier * 0.5 * pi")
		->required();
	dist_ggx_cmd->add_option("-u,--roughness", ggx_params.roughness, "roughness")
		->required();
	dist_ggx_cmd->add_option("-o,--output", ggx_params.out_filename, "name of output data file")
		->default_val("ggx.dat");
	dist_ggx_cmd->callback([&]() {
		std::cout << "resolution = " << ggx_params.resolution << std::endl;
		std::cout << "directory = " << out_dir << std::endl;
		std::cout << "theta multiplier = " << ggx_params.theta_mult << std::endl;
		std::cout << "roughness = " << ggx_params.roughness << std::endl;
		std::cout << "out filename = " << ggx_params.out_filename << std::endl;

		float theta_v = glm::half_pi<float>() * ggx_params.theta_mult;
		glm::vec3 v = UVA::spherical(theta_v, 0.0f);
		float alpha = ggx_params.roughness * ggx_params.roughness;
		std::cout << "alpha = " << alpha << std::endl;
		BRDF brdf(alpha);
		print_ggx(out_dir + "/" + ggx_params.out_filename, brdf, v, ggx_params.resolution);
	});

	// subcommand dist-lct
	auto dist_lct_cmd = app.add_subcommand("dist-lct", "Print values of LCT BRDF on upper hemisphere to a file");
	struct DistLCTParams {
		int resolution = 0;
		std::vector<float> m_params;
		float y_angle;
		float albedo;
		std::string out_filename = "";
	};
	DistLCTParams lct_params;
	dist_lct_cmd->add_option("-r,--res", lct_params.resolution, "sampling resolution of theta and phi")
		->default_val(256);
	dist_lct_cmd->add_option("-m,--m-params", lct_params.m_params, "Matrix m parameters. x-scale, x-z shear, y scale")
		->required()
		->expected(3);
	dist_lct_cmd->add_option("-y,--y-angle", lct_params.y_angle, "Angle of rotation around +y axis, in degrees")
		->default_val(0.0f);
	dist_lct_cmd->add_option("-a,--albedo", lct_params.albedo, "albedo parameter")
		->default_val(1.0f);
	dist_lct_cmd->add_option("-o,--output", lct_params.out_filename, "name of output data file")
		->default_val("lct.dat");
	dist_lct_cmd->callback([&]() {
		std::cout << "resolution = " << lct_params.resolution << std::endl;
		std::cout << "x-scale = " << lct_params.m_params[0] << ", x-z shear = " << lct_params.m_params[1] << ", y scale = " << lct_params.m_params[2] << std::endl;
		std::cout << "y rotation = " << lct_params.y_angle << std::endl;
		std::cout << "albedo = " << lct_params.albedo << std::endl;
		std::cout << "out filename = " << lct_params.out_filename << std::endl;

		glm::quat quat = glm::angleAxis(glm::radians(lct_params.y_angle), glm::vec3(0.0f, 1.0f, 0.0f));
		glm::mat3 rot = glm::toMat3(quat);

		LCT lct(rot,
			{ lct_params.m_params[0],
			lct_params.m_params[1],
			lct_params.m_params[2] },
			lct_params.albedo);
		print_lct(out_dir + "/" + lct_params.out_filename, lct, lct_params.resolution);
	});

	// subcommand: print distribution
	auto dist_cmd = app.add_subcommand("dist", "Print values of LCT and GGX on upper hemisphere to files");
	struct DistParams {
		int resolution = 0;
		float theta_mult = 0.0f;
		float roughness = 0.0f;
		float precond_step = 0.0f;
		std::string lct_filename = "";
		std::string ggx_filename = "";
	};
	DistParams d_params;
	dist_cmd->add_option("-r,--res", d_params.resolution, "sampling resolution of theta and phi")
		->default_val(256);
	dist_cmd->add_option("-t,--theta-multiplier", d_params.theta_mult, "theta multiplier. Viewing angle theta_v = multiplier * 0.5 * pi")
		->required();
	dist_cmd->add_option("-u,--roughness", d_params.roughness, "roughness")
		->required();
	dist_cmd->add_option("-p,--precond-step", d_params.precond_step, "preconditioning step size")
		->default_val(0.015625f); // 1/64
	dist_cmd->add_option("-l,--lct-filename", d_params.lct_filename, "name of LCT data file")
		->default_val("lct.dat");
	dist_cmd->add_option("-g,--ggx-filename", d_params.ggx_filename, "name of GGX BRDF data file")
		->default_val("ggx.dat");
	dist_cmd->callback([&]() {
		std::cout << "resolution = " << d_params.resolution << std::endl;
		std::cout << "directory = " << out_dir << std::endl;
		std::cout << "theta multiplier = " << d_params.theta_mult << std::endl;
		std::cout << "roughness = " << d_params.roughness << std::endl;
		std::cout << "precondition step size = " << d_params.precond_step << std::endl;
		std::cout << "lct filename = " << d_params.lct_filename << std::endl;
		std::cout << "ggx filename = " << d_params.ggx_filename << std::endl;

		// theta_v = 0.99 & alpha = 0.03 works 
		float target_alpha = d_params.roughness * d_params.roughness;
		std::cout << "target alpha = " << target_alpha << std::endl;
		float target_theta_v = glm::half_pi<float>() * d_params.theta_mult;
		glm::vec3 target_v = UVA::spherical(target_theta_v, 0.0f);

		Vector<float, 3> start_m_params = { 1.0f, 0.0f, 1.0f };
		for (float r = 1.0f; r > d_params.roughness; r -= d_params.precond_step) {
			std::cout << "fitting roughness " << r << std::endl;
			float alpha = r * r;
			BRDF brdf(alpha);
			BRDF::AlbedoContext a_ctx = brdf.albedo_2(glm::vec3(0.0f, 0.0f, 1.0f));
			FitResult f_res = fit(glm::vec3(0.0f, 0.0f, 1.0f), brdf, a_ctx.albedo, a_ctx.l_avg, start_m_params);
			start_m_params = f_res.opti_m_params;
		}
		for (float t = 0.0f; t < d_params.theta_mult; t += d_params.precond_step) {
			std::cout << "fitting viewing angle multiplier " << t << std::endl;
			float theta_v = glm::half_pi<float>() * t;
			glm::vec3 v = UVA::spherical(theta_v, 0.0f);
			BRDF brdf(target_alpha);
			BRDF::AlbedoContext a_ctx = brdf.albedo_2(v);
			FitResult f_res = fit(v, brdf, a_ctx.albedo, a_ctx.l_avg, start_m_params);
			start_m_params = f_res.opti_m_params;
		}
		BRDF brdf(target_alpha);
		BRDF::AlbedoContext a_ctx = brdf.albedo_2(target_v);
		FitResult f_res = fit(target_v, brdf, a_ctx.albedo, a_ctx.l_avg, start_m_params);

		LCT lct(f_res.rotation, f_res.opti_m_params, f_res.albedo);

		print_ggx(out_dir + "/" + d_params.ggx_filename, brdf, target_v, d_params.resolution);
		print_lct(out_dir + "/" + d_params.lct_filename, lct, d_params.resolution);
	});

	// subcommand: generate LCT parameter LUTs
	auto gen_luts_cmd = app.add_subcommand("gen_luts", "write matrix M and albedo to textures");
	std::tuple<int, float, float> theta_range_tuple{ 64, 0.0f, 0.99f };
	gen_luts_cmd->add_option("-t,--theta-multiplier",
		theta_range_tuple,
		"range of theta multiplier: <resolution> <min> <max>")
		->default_str("64, 0.0, 0.99");

	std::tuple<int, float, float> roughness_range_tuple{ 64, 0.03f, 1.0f };
	gen_luts_cmd->add_option("-r,--roughness",
		roughness_range_tuple,
		"range of roughness: <resolution> <min> <max>")
		->default_str("64, 0.03, 1.0");

	gen_luts_cmd->callback([&]() {
		TabRange theta_range{ std::get<0>(theta_range_tuple), std::get<1>(theta_range_tuple) , std::get<2>(theta_range_tuple) };
		TabRange roughness_range{ std::get<0>(roughness_range_tuple), std::get<1>(roughness_range_tuple), std::get<2>(roughness_range_tuple) };
		std::vector<LUT> luts = std::move(gen_lut(roughness_range, theta_range));

		// lut_m, lut_b, lut_bm_0, lut_bm_1, lut_inv_0, lut_inv_1, lut_norm_inv

		luts[0].write_to_exr("./out/lut_m.exr");
		luts[1].write_to_exr("./out/lut_b.exr");
		luts[2].write_to_exr("./out/lut_bm_0.exr");
		luts[3].write_to_exr("./out/lut_bm_1.exr");
		luts[4].write_to_exr("./out/lut_inv_0.exr");
		luts[5].write_to_exr("./out/lut_inv_1.exr");
		luts[6].write_to_exr("./out/lut_norm_inv.exr");
	});

	// subcommand: sample LCT parameter LUTs
	auto sample_luts_cmd = app.add_subcommand("sample_luts", "sample textures containing M parameters");
	sample_luts_cmd->add_option("-r,--res", d_params.resolution, "sampling resolution of theta and phi")
		->default_val(256);
	sample_luts_cmd->add_option("-t,--theta-multiplier", d_params.theta_mult, "theta multiplier. Viewing angle theta_v = multiplier * 0.5 * pi")
		->required();
	sample_luts_cmd->add_option("-u,--roughness", d_params.roughness, "roughness")
		->required();
	sample_luts_cmd->add_option("-l,--lct-filename", d_params.lct_filename, "name of LCT data file")
		->default_val("lct.dat");

	std::string filename_0;
	std::string filename_1;
	sample_luts_cmd->add_option("input_0", filename_0, "LUT texture containing first 3 parameters of matrix M")
		->required()
		->check(CLI::ExistingFile);
	sample_luts_cmd->add_option("input_1", filename_1, "LUT texture containing last 2 parameters of matrix M and albedo")
		->required()
		->check(CLI::ExistingFile);

	sample_luts_cmd->callback([&]() {
		LUT lut_bm_0(filename_0);
		LUT lut_bm_1(filename_1);

		float r = d_params.roughness;  // roughness
		float t = d_params.theta_mult; // theta multiplier
		glm::vec4 p0 = lut_bm_0.sample(glm::vec2(r, t));
		glm::vec4 p1 = lut_bm_1.sample(glm::vec2(r, t));
		LCT lct({ p0[0], p0[1], p0[2], p1[0], p1[1] }, p1[2]);

		print_lct(out_dir + "/" + d_params.lct_filename, lct, d_params.resolution);
	});

	app.require_subcommand(1);  // At least 1 subcommand required

	CLI11_PARSE(app, argc, argv);

	return 0;
}
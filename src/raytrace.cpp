#include "scene.h"

#include <thread>
#include "ext/yocto_utils.h"

ray3f eval_camera(const camera* cam, const vec2f& uv) {
    auto h = 2 * std::tan(cam->fovy / 2);
    auto w = h * cam->aspect;
    auto o = cam->frame.o;
    auto q = cam->frame.x * w * (uv.x - 0.5f) +
             cam->frame.y * h * (uv.y - 0.5f) - cam->frame.z + cam->frame.o;
    return ray3f{o, normalize(q - o), ray_eps, 100000};
}


vec4f shade(const scene* scn, const std::vector<instance*>& lights,
    const vec3f& amb, const ray3f& ray) {
    auto isec = intersect_first(scn, ray);
    if (!isec.hit()) return {0, 0, 0, 0};

    auto ist = isec.ist;
    auto shp = ist->shp;
    auto mat = ist->mat;

    auto pos = eval_pos(ist, isec.ei, isec.ew);
    auto norm = eval_norm(ist, isec.ei, isec.ew);
    auto texcoord = eval_texcoord(shp, isec.ei, isec.ew);

    auto ke = mat->ke * eval_texture(mat->ke_txt, texcoord, true);
    auto kd = mat->kd * eval_texture(mat->kd_txt, texcoord, true);
    auto ks = mat->ks * eval_texture(mat->ks_txt, texcoord, true);
    auto kr = mat->kr * eval_texture(mat->kr_txt, texcoord, true);
    auto ns = (mat->rs) ? 2 / std::pow(mat->rs, 4.0f) - 2 : 1e6f;
    auto op = mat->op;

    auto l = ke + kd * amb;
    for (auto lgt : lights) {
        auto lshp = lgt->shp;
        for (auto p : lgt->shp->points) {
            auto lpos = transform_point(lgt->frame, lshp->pos[p]);
            auto lr = length(lpos - pos);
            auto wi = normalize(lpos - pos);
#if 0
            if (intersect_any(scn, ray3f{pos, wi, 0.01f, lr - 0.04f})) continue;
#else
            auto spos = pos;
            auto tr = 1.0f;
            auto sisec = intersect_first(
                scn, {spos, wi, 0.01f, length(lpos - spos) - 0.04f});
            while (sisec.hit() && tr > 0) {
                tr *= 1 - sisec.ist->mat->op;
                spos = eval_pos(sisec.ist, sisec.ei, sisec.ew);
                sisec = intersect_first(
                    scn, {spos, wi, 0.01f, length(lpos - spos) - 0.04f});
            }
#endif

            auto le = lgt->mat->ke * tr / (lr * lr);
            auto wo = -ray.d;
            auto wh = normalize(wi + wo);
            if (!shp->triangles.empty()) {
                l += le * kd * max(0.0f, dot(norm, wi)) +
                     le * ks * pow(max(0.0f, dot(norm, wh)), ns);
            } else if (!shp->lines.empty()) {
                l += le * kd *
                         sqrt(clamp(
                             1 - dot(norm, wi) * dot(norm, wi), 0.0f, 1.0f)) +
                     le * ks *
                         pow(sqrt(clamp(1 - dot(norm, wh) * dot(norm, wh), 0.0f,
                                 1.0f)),
                             ns);
            }
        }
    }

    if (kr != vec3f{0, 0, 0}) {
        auto rray = ray3f{pos, ray.d - norm * 2 * dot(ray.d, norm), 0.001f};
        auto cr = shade(scn, lights, amb, rray);
        l += kr * vec3f{cr.x, cr.y, cr.z};
    }

    if (op < 1) {
        auto rray = ray3f{pos, ray.d, 0.001f};
        auto cr = shade(scn, lights, amb, rray);
        l = l * op + vec3f{cr.x, cr.y, cr.z} * (1 - op);
    }

    return {l.x, l.y, l.z, 1};
}


void
raytrace_region(const scene* scn, vec3f amb, int samples,
                const camera* cam, image4f& img, const std::vector<instance*>& lights,
                int x0, int y0, int x1, int y1)
{
    for (auto j = y0; j < y1; j++)
    {
        for (auto i = x0; i < x1; i++)
        {
            img.at(i, j) = {0, 0, 0, 0};
            for (auto sj = 0; sj < samples; sj++)
            {
                for (auto si = 0; si < samples; si++)
                {
                    auto u = (i + (si + 0.5f) / samples) / img.width;
                    auto v = ((img.height - j) + (sj + 0.5f) / samples) / img.height;
                    auto ray = eval_camera(cam, {u, v});
                    img.at(i, j) += shade(scn, lights, amb, ray);
                }
            }
            img.at(i, j) *= 1.0f / (float)(samples * samples);
        }
    }
}


image4f raytrace(const scene* scn, const vec3f& amb, int resolution, int samples)
{
    auto cam = scn->cameras.front();
    auto img = image4f((int)std::round(cam->aspect * resolution), resolution);

    auto lights = std::vector<instance*>();
    for (auto ist : scn->instances) {
        if (ist->mat->ke == vec3f{0, 0, 0}) continue;
        if (ist->shp->points.empty()) continue;
        lights.push_back(ist);
    }
#if 0
    for (auto j = 0; j < img.height; j++) {
        for (auto i = 0; i < img.width; i++) {
            img.at(i, j) = {0, 0, 0, 0};
            for (auto sj = 0; sj < samples; sj++) {
                for (auto si = 0; si < samples; si++) {
                    auto u = (i + (si + 0.5f) / samples) / img.width;
                    auto v =
                        ((img.height - j) + (sj + 0.5f) / samples) / img.height;
                    auto ray = eval_camera(cam, {u, v});
                    img.at(i, j) += shade(scn, lights, amb, ray);
                }
            }
            img.at(i, j) *= 1.0f / (float)(samples * samples);
        }
    }
#else
    std::thread* threads       = NULL;
    unsigned int threads_count = std::thread::hardware_concurrency();
    if(threads_count > 1)
        threads = new std::thread[threads_count-1];
    
    int height_per_thread = img.height / threads_count;
    int leftover_height   = img.height % threads_count;
    
    for(unsigned int ti = 1; ti < threads_count; ++ti)
    {
        int y0 = ti * height_per_thread;
        int y1 = y0 + height_per_thread;
        if(leftover_height > 0)
            ++y1, --leftover_height;
        
        threads[ti-1] = std::thread(raytrace_region,
                                    scn, amb, samples,
                                    cam, std::ref(img), std::ref(lights), 
                                    0, y0, img.width, y1);
    }
    
    raytrace_region(scn, amb, samples,
                    cam, img, lights,
                    0, 0, img.width, height_per_thread);
    

    for(unsigned int ti = 1; ti < threads_count; ++ti)
    {
        threads[ti-1].join();
    }
    
    if(threads)
        delete[] threads;
#endif

    return img;
}


int main(int argc, char** argv) {
    // command line parsing
    auto parser =
        yu::cmdline::make_parser(argc, argv, "raytrace", "raytrace scene");
    auto resolution = yu::cmdline::parse_opti(
        parser, "--resolution", "-r", "vertical resolution", 720);
    auto samples = yu::cmdline::parse_opti(
        parser, "--samples", "-s", "per-pixel samples", 1);
    auto amb = yu::cmdline::parse_optf(
        parser, "--ambient", "-a", "ambient color", 0.1f);
    auto imageout = yu::cmdline::parse_opts(
        parser, "--output", "-o", "output image", "out.png");
    auto scenein = yu::cmdline::parse_args(
        parser, "scenein", "input scene", "scene.obj", true);
    yu::cmdline::check_parser(parser);

    // load scene
    printf("loading scene %s\n", scenein.c_str());
    auto scn = load_scene(scenein);

    // create bvh
    printf("creating bvh\n");
    build_bvh(scn, false);

    // raytrace
    printf("tracing scene\n");
    auto hdr = raytrace(scn, vec3f{amb, amb, amb}, resolution, samples);

    // tonemap and save
    printf("saving image %s\n", imageout.c_str());
    save_hdr_or_ldr(imageout, hdr);
}

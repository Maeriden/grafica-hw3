
#include "scene.h"

#ifndef YSCN_NO_IMAGE
    #define YSCN_NO_IMAGE
#endif
#include "ext/yocto_scn.h"

#include <map>
#include <set>

// LITC begin smooth_normals
void compute_smooth_normals(shape* shp) {
    shp->norm.resize(shp->pos.size());
    for (auto& norm : shp->norm) norm = {0, 0, 0};

    for (auto l : shp->lines) {
        auto n = line_tangent(shp->pos[l.x], shp->pos[l.y]);
        auto w = line_length(shp->pos[l.x], shp->pos[l.y]);
        shp->norm[l.x] += n * w;
        shp->norm[l.y] += n * w;
    }

    for (auto t : shp->triangles) {
        auto n = triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        auto w = triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        shp->norm[t.x] += n * w;
        shp->norm[t.y] += n * w;
        shp->norm[t.z] += n * w;
    }

    for (auto t : shp->quads) {
        auto n = triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]) +
                 triangle_normal(shp->pos[t.x], shp->pos[t.z], shp->pos[t.w]);
        auto w = triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]) +
                 triangle_area(shp->pos[t.x], shp->pos[t.z], shp->pos[t.w]);
        shp->norm[t.x] += n * w;
        shp->norm[t.y] += n * w;
        shp->norm[t.z] += n * w;
        shp->norm[t.w] += n * w;
    }

    for (auto& norm : shp->norm) norm = normalize(norm);
}
// LITC end smooth_normals

// LITC begin animfuncs
void update_animation(instance* ist, float time) {
    auto anim = ist->anim;
    if (!anim) return;
    auto idx = int(time / anim->delta_t) % anim->frame_keyframes.size();
    ist->frame = anim->frame_keyframes[idx];
    if (!anim->pos_keyframes.empty()) ist->shp->pos = anim->pos_keyframes[idx];
    if (!anim->norm_keyframes.empty())
        ist->shp->norm = anim->norm_keyframes[idx];
}

void add_keyframe(instance* ist) {
    ist->anim->frame_keyframes.push_back(ist->frame);
    ist->anim->pos_keyframes.push_back(ist->shp->pos);
    ist->anim->norm_keyframes.push_back(ist->shp->norm);
}
// LITC end animfuncs

// triangle tangent and bitangent from uv (not othornormalized with themselfves
// not the normal)
// see http://www.terathon.com/code/tangent.html and
// https://gist.github.com/aras-p/2843984
// LITC begin compute_tangent_space
std::pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
    // normal points up from texture space
    auto p = v1 - v0;
    auto q = v2 - v0;
    auto s = vec2f{uv1.x - uv0.x, uv2.x - uv0.x};
    auto t = vec2f{uv1.y - uv0.y, uv2.y - uv0.y};
    auto div = s.x * t.y - s.y * t.x;

    if (div > 0) {
        auto tu = vec3f{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                      t.y * p.z - t.x * q.z} /
                  div;
        auto tv = vec3f{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                      s.x * q.z - s.y * p.z} /
                  div;
        return {tu, tv};
    } else {
        return {{1, 0, 0}, {0, 1, 0}};
    }
}

void compute_tangent_space(shape* shp) {
    auto tangu = std::vector<vec3f>(shp->pos.size(), {0, 0, 0});
    auto tangv = std::vector<vec3f>(shp->pos.size(), {0, 0, 0});

    for (auto t : shp->triangles) {
        auto tutv = triangle_tangents_fromuv(shp->pos[t.x], shp->pos[t.y],
            shp->pos[t.z], shp->texcoord[t.x], shp->texcoord[t.y],
            shp->texcoord[t.z]);
        auto w = triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        tangu[t.x] += tutv.first * w;
        tangu[t.y] += tutv.first * w;
        tangu[t.z] += tutv.first * w;
        tangv[t.x] += tutv.second * w;
        tangv[t.y] += tutv.second * w;
        tangv[t.z] += tutv.second * w;
    }

    shp->tangsp.resize(shp->pos.size());
    for (auto i = 0; i < tangu.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], shp->norm[i]);
        auto s =
            (dot(cross(shp->norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        shp->tangsp[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
}
// LITC end compute_tangent_space

// LITC begin load
inline frame3f to_frame(const ym::frame3f& m) {
    return {{m.x.x, m.x.y, m.x.z}, {m.y.x, m.y.y, m.y.z}, {m.z.x, m.z.y, m.z.z},
        {m.o.x, m.o.y, m.o.z}};
}

scene* load_scene(const std::string& filename) {
    auto yscn_ = (yscn::scene*)nullptr;
    try {
        auto load_opts = yscn::load_options();
        load_opts.load_textures = false;
        yscn_ = yscn::load_scene(filename, load_opts);
    } catch (std::exception& e) {
        printf("could not load scene\n");
        exit(1);
    }

    auto add_opts = yscn::add_elements_options::none();
    add_opts.shape_instances = true;
    add_opts.default_names = true;
    add_opts.default_camera = true;
    add_opts.pointline_radius = 0.001f;
    add_opts.texture_data = false;
    yscn::add_elements(yscn_, add_opts);

    auto scn = new scene();

    for (auto ycam : yscn_->cameras) {
        auto cam = new camera();
        cam->name = ycam->name;
        cam->frame = to_frame(ycam->frame);
        cam->fovy = ycam->yfov;
        cam->aspect = ycam->aspect;
        cam->focus = ycam->focus;
        cam->aperture = ycam->aperture;
        scn->cameras.push_back(cam);
    }

    auto dirname = std::string();
    if (filename.rfind("/") != filename.npos) {
        dirname = filename.substr(0, filename.rfind("/") + 1);
    }
    auto textures = std::map<yscn::texture*, texture*>{{nullptr, nullptr}};
    for (auto ytxt : yscn_->textures) {
        auto txt = new texture();
        txt->filename = ytxt->path;
        if (txt->filename.substr(txt->filename.length() - 4) == ".hdr") {
            txt->hdr = load_image4f(dirname + ytxt->path);
        } else {
            txt->ldr = load_image4b(dirname + ytxt->path);
        }
        scn->textures.push_back(txt);
        textures[ytxt] = txt;
    }

    auto materials = std::map<yscn::material*, material*>();
    for (auto ymat : yscn_->materials) {
        auto mat = new material();
        mat->name = ymat->name;
        mat->ke = {ymat->ke[0], ymat->ke[1], ymat->ke[2]};
        mat->kd = {ymat->kd[0], ymat->kd[1], ymat->kd[2]};
        mat->ks = {ymat->ks[0], ymat->ks[1], ymat->ks[2]};
        mat->rs = ymat->rs;
        mat->kr = {ymat->kr[0], ymat->kr[1], ymat->kr[2]};
        mat->op = ymat->op;
        mat->ke_txt = textures.at(ymat->ke_txt.txt);
        mat->kd_txt = textures.at(ymat->kd_txt.txt);
        mat->ks_txt = textures.at(ymat->ks_txt.txt);
        mat->kr_txt = textures.at(ymat->kr_txt.txt);
        mat->norm_txt = textures.at(ymat->norm_txt.txt);
        mat->disp_txt = textures.at(ymat->disp_txt.txt);
        mat->disp_height = 0.01f;
        scn->materials.push_back(mat);
        materials[ymat] = mat;
    }

    auto shapes = std::map<yscn::shape*, shape*>();
    for (auto yshp : yscn_->shapes) {
        auto shp = new shape();
        shp->name = yshp->name;
        shp->pos = *(std::vector<vec3f>*)&yshp->pos;
        shp->norm = *(std::vector<vec3f>*)&yshp->norm;
        shp->texcoord = *(std::vector<vec2f>*)&yshp->texcoord;
        shp->radius = yshp->radius;
        shp->points = yshp->points;
        shp->lines = *(std::vector<vec2i>*)&yshp->lines;
        shp->triangles = *(std::vector<vec3i>*)&yshp->triangles;
        scn->shapes.push_back(shp);
        shapes[yshp] = shp;
    }

    for (auto yist : yscn_->instances) {
        auto yshp = yist->shp;
        auto ist = new instance();
        ist->name = yist->name;
        ist->frame = to_frame(yist->frame);
        ist->mat = materials[yshp->mat];
        ist->shp = shapes[yshp];
        scn->instances.push_back(ist);
    }

    for (auto yenv : yscn_->environments) {
        auto ymat = yenv->mat;
        auto env = new environment();
        env->name = yenv->name;
        env->frame = to_frame(yenv->frame);
        env->ke = {ymat->ke[0], ymat->ke[1], ymat->ke[2]};
        env->ke_txt = textures.at(ymat->ke_txt.txt);
        scn->environments.push_back(env);
    }

    for (auto ist : scn->instances) {
        auto shp = ist->shp;
        if (shp->norm.empty()) compute_smooth_normals(shp);
        if (shp->tangsp.empty() && ist->mat->norm_txt)
            compute_tangent_space(shp);
    }

    return scn;
}
// LITC end load

// LITC begin save
inline ym::frame3f to_frame(const frame3f& m) {
    return {{m.x.x, m.x.y, m.x.z}, {m.y.x, m.y.y, m.y.z}, {m.z.x, m.z.y, m.z.z},
        {m.o.x, m.o.y, m.o.z}};
}

void save_scene(const std::string& filename, const scene* scn) {
    auto yscn_ = new yscn::scene();

    for (auto cam : scn->cameras) {
        auto ycam = new yscn::camera();
        ycam->name = cam->name;
        ycam->frame = to_frame(cam->frame);
        ycam->yfov = cam->fovy;
        ycam->aspect = cam->aspect;
        ycam->focus = cam->focus;
        ycam->aperture = cam->aperture;
        yscn_->cameras.push_back(ycam);
    }

    auto dirname = std::string();
    if (filename.rfind("/") != filename.npos) {
        dirname = filename.substr(0, filename.rfind("/") + 1);
    }
    auto textures = std::map<texture*, yscn::texture*>{{nullptr, nullptr}};
    for (auto txt : scn->textures) {
        auto ytxt = new yscn::texture();
        ytxt->path = txt->filename;
        yscn_->textures.push_back(ytxt);
        textures[txt] = ytxt;
    }

    auto materials = std::map<material*, yscn::material*>();
    for (auto mat : scn->materials) {
        auto ymat = new yscn::material();
        ymat->name = mat->name;
        ymat->ke = {mat->ke.x, mat->ke.y, mat->ke.z};
        ymat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
        ymat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
        ymat->rs = mat->rs;
        ymat->kr = {mat->kr.x, mat->kr.y, mat->kr.z};
        ymat->op = mat->op;
        ymat->ke_txt.txt = textures.at(mat->ke_txt);
        ymat->kd_txt.txt = textures.at(mat->kd_txt);
        ymat->ks_txt.txt = textures.at(mat->ks_txt);
        ymat->kr_txt.txt = textures.at(mat->kr_txt);
        ymat->norm_txt.txt = textures.at(mat->norm_txt);
        ymat->disp_txt.txt = textures.at(mat->disp_txt);
        ymat->disp_txt.scale = 0.01f;
        yscn_->materials.push_back(ymat);
        materials[mat] = ymat;
    }

    auto shapes = std::map<shape*, yscn::shape*>();
    for (auto shp : scn->shapes) {
        auto yshp = new yscn::shape();
        yshp->name = shp->name;
        yshp->pos = *(std::vector<ym::vec3f>*)&shp->pos;
        yshp->norm = *(std::vector<ym::vec3f>*)&shp->norm;
        yshp->texcoord = *(std::vector<ym::vec2f>*)&shp->texcoord;
        yshp->radius = shp->radius;
        yshp->points = shp->points;
        yshp->lines = *(std::vector<ym::vec2i>*)&shp->lines;
        yshp->triangles = *(std::vector<ym::vec3i>*)&shp->triangles;
        yshp->quads = *(std::vector<ym::vec4i>*)&shp->quads;
        yscn_->shapes.push_back(yshp);
        shapes[shp] = yshp;
    }

    for (auto ist : scn->instances) {
        auto yist = new yscn::instance();
        yist->name = ist->name;
        yist->frame = to_frame(ist->frame);
        yist->shp = shapes[ist->shp];
        yist->shp->mat = materials[ist->mat];
        yscn_->instances.push_back(yist);
    }

    for (auto env : scn->environments) {
        auto yenv = new yscn::environment();
        yenv->name = env->name;
        yenv->frame = to_frame(env->frame);
        yenv->mat = new yscn::material();
        yenv->mat->ke = {env->ke.x, env->ke.y, env->ke.z};
        yenv->mat->ke_txt.txt = textures.at(env->ke_txt);
        yscn_->environments.push_back(yenv);
    }

    try {
        auto save_opts = yscn::save_options();
        save_opts.save_textures = false;
        yscn::save_scene(filename, yscn_, save_opts);
    } catch (std::exception& e) {
        printf("could not load scene\n");
        exit(1);
    }
}
// LITC end save

// LITC begin intersect_triangle
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& dist, vec4f& ew) {
    // compute triangle edges
    auto e1 = v1 - v0;
    auto e2 = v2 - v0;

    // compute denominator
    auto r = cross(ray.d, e2);
    auto den = dot(r, e1);

    // check denominator and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (den == 0) return false;
    float inv_den = 1.0f / den;

    // compute and check first barycentric coordinate
    auto c = ray.o - v0;
    auto w1 = dot(r, c) * inv_den;
    if (w1 < 0 || w1 > 1) return false;

    // compute and check second barycentric coordinate
    auto s = cross(c, e1);
    auto w2 = dot(s, ray.d) * inv_den;
    if (w2 < 0.0 || w1 + w2 > 1.0) return false;

    // compute and check ray parameter
    auto t = dot(s, e2) * inv_den;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    dist = t;
    ew = {1 - w1 - w2, w1, w2, 0};

    return true;
}
// LITC end intersect_triangle

// LITC begin intersect_point
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& dist, vec4f& ew) {
    auto w = p - ray.o;
    auto t = dot(w, ray.d) / dot(ray.d, ray.d);
    if (t < ray.tmin || t > ray.tmax) return false;

    auto rp = eval_ray(ray, t);
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    dist = t;
    ew = {1, 0, 0, 0};

    return true;
}
// LITC end intersect_point

// LITC begin intersect_line
inline bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& dist, vec4f& ew) {
    auto u = ray.d, v = v1 - v0, w = ray.o - v0;

    auto a = dot(u, u), b = dot(u, v), c = dot(v, v), d = dot(u, w),
         e = dot(v, w);
    auto det = a * c - b * b;
    if (det == 0) return false;

    auto t = (b * e - c * d) / det, s = (a * e - b * d) / det;
    if (t < ray.tmin || t > ray.tmax) return false;

    s = clamp(s, (float)0, (float)1);
    auto p0 = eval_ray(ray, t), p1 = eval_ray(ray3f{v0, v1 - v0}, s);
    auto p01 = p0 - p1;
    auto r = r0 * (1 - s) + r1 * s;
    if (dot(p01, p01) > r * r) return false;

    dist = t;
    ew = {1 - s, s, 0, 0};

    return true;
}
// LITC end intersect_line

// LITC begin intersect_ray_shape
bool intersect_ray(
    const shape* shp, const ray3f& ray, float& dist, int& iei, vec4f& iew) {
    auto hit = false;
    auto tray = ray;
    for (auto ei = 0; ei < shp->triangles.size(); ei++) {
        auto e = shp->triangles[ei];
        if (!intersect_triangle(
                tray, shp->pos[e.x], shp->pos[e.y], shp->pos[e.z], dist, iew))
            continue;
        iei = ei;
        tray.tmax = dist;
        hit = true;
    }
    for (auto ei = 0; ei < shp->lines.size(); ei++) {
        auto e = shp->lines[ei];
        if (!intersect_line(tray, shp->pos[e.x], shp->pos[e.y],
                shp->radius[e.x], shp->radius[e.y], dist, iew))
            continue;
        iei = ei;
        tray.tmax = dist;
        hit = true;
    }
    for (auto ei = 0; ei < shp->points.size(); ei++) {
        auto e = shp->points[ei];
        if (!intersect_point(tray, shp->pos[e], shp->radius[e], dist, iew))
            continue;
        iei = ei;
        tray.tmax = dist;
        hit = true;
    }
    return hit;
}
// LITC end intersect_ray_shape

// LITC begin intersect_ray_scene
bool intersect_ray(const scene* scn, const ray3f& ray, float& dist,
    instance*& iist, int& iei, vec4f& iew) {
    auto hit = false;
    auto tray = ray;
    for (auto ist : scn->instances) {
        auto lray = transform_ray_inverse(ist->frame, tray);
        if (!intersect_ray(ist->shp, lray, dist, iei, iew)) continue;
        tray.tmax = dist;
        iist = ist;
        hit = true;
    }
    return hit;
}
// LITC end intersect_ray_scene

// LITC begin intersect_ray_easy
intersection3f intersect_scene(const scene* scn, const ray3f& ray) {
    intersection3f isec;
    if (!intersect_ray(scn, ray, isec.dist, isec.ist, isec.ei, isec.ew))
        return {};
    return isec;
}
// LITC end intersect_ray_easy

// LITC begin intersect_bbox
inline bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox) {
    auto invd = vec3f{1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z};
    auto t0 = (bbox.min - ray.o) * invd;
    auto t1 = (bbox.max - ray.o) * invd;
    if (invd.x < 0) std::swap(t0.x, t1.x);
    if (invd.y < 0) std::swap(t0.y, t1.y);
    if (invd.z < 0) std::swap(t0.z, t1.z);
    auto tmin = max(t0.z, max(t0.y, max(t0.x, ray.tmin)));
    auto tmax = min(t1.z, min(t1.y, min(t1.x, ray.tmax)));
    tmax *= 1.00000024f;
    return tmin <= tmax;
}
// LITC end intersect_bbox

// LITC begin intersect_bvh_shape
bool intersect_bvh(const shape* shp, const ray3f& ray, bool any, float& dist,
    int& ei, vec4f& ew) {
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    auto bvh = shp->bvh;

    auto tray = ray;
    auto hit = false;

    while (node_cur) {
        auto node = bvh->nodes[node_stack[--node_cur]];
        if (!intersect_check_bbox(tray, node.bbox)) continue;

        if (!node.isleaf) {
            for (auto i = node.start; i < node.start + node.count; i++) {
                node_stack[node_cur++] = i;
            }
        } else if (!shp->triangles.empty()) {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto e = shp->triangles[bvh->leaf_prims[i]];
                if (!intersect_triangle(tray, shp->pos[e.x], shp->pos[e.y],
                        shp->pos[e.z], dist, ew))
                    continue;
                hit = true;
                tray.tmax = dist;
                ei = bvh->leaf_prims[i];
                if (any) return true;
            }
        } else if (!shp->lines.empty()) {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto e = shp->lines[bvh->leaf_prims[i]];
                if (!intersect_line(tray, shp->pos[e.x], shp->pos[e.y],
                        shp->radius[e.x], shp->radius[e.y], dist, ew))
                    continue;
                hit = true;
                tray.tmax = dist;
                ei = bvh->leaf_prims[i];
                if (any) return true;
            }
        } else if (!shp->points.empty()) {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto e = shp->points[bvh->leaf_prims[i]];
                if (!intersect_point(
                        tray, shp->pos[e], shp->radius[e], dist, ew))
                    continue;
                hit = true;
                tray.tmax = dist;
                ei = bvh->leaf_prims[i];
                if (any) return true;
            }
        }
    }

    return hit;
}
// LITC end intersect_bvh_shape

// LITC begin intersect_bvh_scene
bool intersect_bvh(const scene* scn, const ray3f& ray, bool any, float& dist,
    instance*& ist, int& ei, vec4f& ew) {
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    auto bvh = scn->bvh;

    auto tray = ray;
    auto hit = false;

    while (node_cur) {
        auto node = bvh->nodes[node_stack[--node_cur]];
        if (!intersect_check_bbox(tray, node.bbox)) continue;

        if (!node.isleaf) {
            for (auto i = node.start; i < node.start + node.count; i++) {
                node_stack[node_cur++] = i;
            }
        } else {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto is = scn->instances[bvh->leaf_prims[i]];
                auto lray = transform_ray_inverse(is->frame, tray);
                if (!intersect_bvh(is->shp, lray, any, dist, ei, ew)) continue;
                tray.tmax = dist;
                ist = is;
                hit = true;
                if (any) return hit;
            }
        }
    }

    return hit;
}
// LITC end intersect_bvh_scene

// LITC begin intersect_bvh_easy
intersection3f intersect_first(const scene* scn, const ray3f& ray) {
    intersection3f isec;
    if (!intersect_bvh(scn, ray, false, isec.dist, isec.ist, isec.ei, isec.ew))
        return {};
    return isec;
}
bool intersect_any(const scene* scn, const ray3f& ray) {
    intersection3f isec;
    return (bool)intersect_bvh(
        scn, ray, true, isec.dist, isec.ist, isec.ei, isec.ew);
}
// LITC end intersect_bvh_easy

// LITC begin intersect_bvh_other_easy
intersection3f intersect_first(const shape* shp, const ray3f& ray) {
    intersection3f isec;
    if (!intersect_bvh(shp, ray, false, isec.dist, isec.ei, isec.ew)) return {};
    return isec;
}
bool intersect_any(const shape* shp, const ray3f& ray) {
    intersection3f isec;
    return (bool)intersect_bvh(shp, ray, true, isec.dist, isec.ei, isec.ew);
}
// LITC end intersect_bvh_other_easy

// LITC begin bound_prim
struct bound_prim {
    bbox3f bbox;   // bounding box
    vec3f center;  // center
    int pid;       // primitive id
};
// LITC end bound_prim

// LITC begin build_bvh_decl
bvh_tree* build_bvh(std::vector<bound_prim>& bound_prims, bool equal_num);
// LITC end build_bvh_decl

// LITC begin build_bvh_shape
bbox3f expand_bbox(const bbox3f& bbox, const vec3f& p, float r) {
    return expand_bbox(bbox, {p - vec3f{r, r, r}, p + vec3f{r, r, r}});
}

void build_bvh(shape* shp, bool equal_num) {
    auto bound_prims = std::vector<bound_prim>();
    for (auto ei = 0; ei < shp->points.size(); ei++) {
        auto e = shp->points[ei];
        auto bbox = invalid_bbox3f;
        bbox = expand_bbox(bbox, shp->pos[e], shp->radius[e]);
        bound_prims.push_back({bbox, (bbox.min + bbox.max) / 2.0f, ei});
    }
    for (auto ei = 0; ei < shp->lines.size(); ei++) {
        auto e = shp->lines[ei];
        auto bbox = invalid_bbox3f;
        bbox = expand_bbox(bbox, shp->pos[e.x], shp->radius[e.x]);
        bbox = expand_bbox(bbox, shp->pos[e.y], shp->radius[e.y]);
        bound_prims.push_back({bbox, (bbox.min + bbox.max) / 2.0f, ei});
    }
    for (auto ei = 0; ei < shp->triangles.size(); ei++) {
        auto e = shp->triangles[ei];
        auto bbox = invalid_bbox3f;
        bbox = expand_bbox(bbox, shp->pos[e.x], 0);
        bbox = expand_bbox(bbox, shp->pos[e.y], 0);
        bbox = expand_bbox(bbox, shp->pos[e.z], 0);
        bound_prims.push_back({bbox, (bbox.min + bbox.max) / 2.0f, ei});
    }

    shp->bvh = build_bvh(bound_prims, equal_num);
}
// LITC end build_bvh_shape

// LITC begin build_bvh_scene
void build_bvh(scene* scn, bool equal_num) {
    for (auto shp : scn->shapes) { build_bvh(shp, equal_num); }

    auto bound_prims = std::vector<bound_prim>();
    for (auto ii = 0; ii < scn->instances.size(); ii++) {
        auto ist = scn->instances[ii];
        auto bbox = bbox_to_world(ist->frame, ist->shp->bvh->nodes[0].bbox);
        bound_prims.push_back({bbox, (bbox.min + bbox.max) / 2.0f, ii});
    }

    scn->bvh = build_bvh(bound_prims, equal_num);
}
// LITC end build_bvh_scene

// LITC begin make_node
bool split_prims(std::vector<bound_prim>& sorted_prim, int start, int end,
    bool equalnum, int& axis, int& mid);

void make_node(bvh_tree* bvh, int nid, std::vector<bound_prim>& leaf_prims,
    int start, int end, bool equal_num) {
    auto node = &bvh->nodes[nid];

    node->bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) {
        node->bbox = expand_bbox(node->bbox, leaf_prims[i].bbox);
    }

    auto split = false;
    auto axis = -1, mid = -1;
    if (end - start > 4) {
        split = split_prims(leaf_prims, start, end, equal_num, axis, mid);
    }

    if (!split) {
        node->isleaf = true;
        node->start = start;
        node->count = end - start;
    } else {
        assert(mid > start && mid < end);
        node->isleaf = false;
        node->axis = axis;
        auto first = (int)bvh->nodes.size();
        node->start = first;
        node->count = 2;
        bvh->nodes.push_back({});
        bvh->nodes.push_back({});
        make_node(bvh, first, leaf_prims, start, mid, equal_num);
        make_node(bvh, first + 1, leaf_prims, mid, end, equal_num);
    }
}
// LITC end make_node

// LITC begin split_prims
bool split_prims(std::vector<bound_prim>& sorted_prim, int start, int end,
    bool equalnum, int& axis, int& mid) {
    auto centroid_bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++)
        centroid_bbox = expand_bbox(centroid_bbox, sorted_prim[i].center);
    auto size = centroid_bbox.max - centroid_bbox.min;

    if (size == vec3f{0, 0, 0}) return false;

    if (size.x >= size.y && size.x >= size.z)
        axis = 0;
    else if (size.y >= size.x && size.y >= size.z)
        axis = 1;
    else
        axis = 2;

    if (equalnum) {
        mid = (start + end) / 2;
        std::nth_element(sorted_prim.begin() + start, sorted_prim.begin() + mid,
            sorted_prim.begin() + end, [axis](auto a, auto b) {
                return (&a.center.x)[axis] < (&b.center.x)[axis];
            });
    } else {
        auto half = (centroid_bbox.min + centroid_bbox.max) / 2;
        mid = (int)(std::partition(sorted_prim.begin() + start,
                        sorted_prim.begin() + end,
                        [axis, half](auto a) {
                            return (&a.center.x)[axis] < (&half.x)[axis];
                        }) -
                    sorted_prim.begin());
    }
    return true;
}
// LITC end split_prims

// LITC begin build_bvh
bvh_tree* build_bvh(std::vector<bound_prim>& bound_prims, bool equal_num) {
    auto bvh = new bvh_tree();

    bvh->nodes.reserve(bound_prims.size() * 2);
    bvh->nodes.push_back({});
    make_node(bvh, 0, bound_prims, 0, bound_prims.size(), equal_num);
    bvh->nodes.shrink_to_fit();

    bvh->leaf_prims.resize(bound_prims.size());
    for (int i = 0; i < bound_prims.size(); i++) {
        bvh->leaf_prims[i] = bound_prims[i].pid;
    }

    return bvh;
}
// LITC end build_bvh

// Tesselates a shape
shape* tesselate(shape* shp) {
    if (!shp->quads.empty()) {
        auto _aux = std::vector<ym::vec4f>();
        ym::tesselate_quads((std::vector<ym::vec4i>&)shp->quads,
            (std::vector<ym::vec3f>&)shp->pos,
            (std::vector<ym::vec3f>&)shp->norm,
            (std::vector<ym::vec2f>&)shp->texcoord, _aux, shp->radius, true);
    }
    if (!shp->triangles.empty()) {
        auto _aux = std::vector<ym::vec4f>();
        ym::tesselate_triangles((std::vector<ym::vec3i>&)shp->triangles,
            (std::vector<ym::vec3f>&)shp->pos,
            (std::vector<ym::vec3f>&)shp->norm,
            (std::vector<ym::vec2f>&)shp->texcoord, _aux, shp->radius, true);
    }
    return shp;
}

vec3f lookup_texture(const texture* txt, int i, int j, bool srgb) {
    if (txt->hdr.width) {
        auto v = txt->hdr.at(i, j);
        return {v.x, v.y, v.z};
    }
    if (txt->ldr.width) {
        auto v = txt->ldr.at(i, j);
        if (srgb)
            return {powf(v.x / 255.f, 2.2f), powf(v.y / 255.f, 2.2f),
                powf(v.z / 255.f, 2.2f)};
        else
            return {v.x / 255.f, v.y / 255.f, v.z / 255.f};
    }
    return {0, 0, 0};
}

vec3f eval_texture(const texture* txt, const vec2f& texcoord, bool srgb) {
    if (!txt) return {1, 1, 1};
    auto width = (txt->hdr.width) ? txt->hdr.width : txt->ldr.width;
    auto height = (txt->hdr.height) ? txt->hdr.height : txt->ldr.height;
    // get coordinates normalized for tiling
    auto s = fmodf(texcoord.x, 1.0f) * width,
         t = fmodf(texcoord.y, 1.0f) * height;
    if (s < 0) s += width;   // since in C++ mod is negative
    if (t < 0) t += height;  // for negative input

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, width), j = clamp((int)t, 0, height);
    auto wi = s - i, wj = t - j;
    auto i1 = (i + 1) % width, j1 = (j + 1) % height;

    // handle interpolation
    return lookup_texture(txt, i, j, srgb) * (1 - wi) * (1 - wj) +
           lookup_texture(txt, i, j1, srgb) * (1 - wi) * wj +
           lookup_texture(txt, i1, j, srgb) * wi * (1 - wj) +
           lookup_texture(txt, i1, j1, srgb) * wi * wj;
}

texture* make_grid_texture(
    int width, int height, int tile, const vec4b& c0, const vec4b& c1) {
    image4b pixels(width, height);
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            pixels.at(i, j) = (c) ? c0 : c1;
        }
    }
    return new texture{"grid", {}, pixels};
}

texture* make_bumpdimple_texture(int width, int height, int tile) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            auto g = (unsigned char)clamp((int)(h * 255), 0, 255);
            pixels.at(i, j) = vec4b{g, g, g, 255};
        }
    }
    return new texture{"bumpdimple", {}, pixels};
}

shape* make_cube(const std::string& name) {
    static auto cube_pos = std::vector<vec3f>{{-1, -1, -1}, {-1, +1, -1},
        {+1, +1, -1}, {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1},
        {+1, -1, +1}};
    static auto cube_quads = std::vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_quad_uv =
        std::vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    auto shp = new shape{name};
    shp->pos = cube_pos;
    shp->quads = cube_quads;
    return shp;
}

shape* make_monkey(const std::string& name) {
    static auto suzanne_pos = std::vector<vec3f>{{0.4375, 0.1640625, 0.765625},
        {-0.4375, 0.1640625, 0.765625}, {0.5, 0.09375, 0.6875},
        {-0.5, 0.09375, 0.6875}, {0.546875, 0.0546875, 0.578125},
        {-0.546875, 0.0546875, 0.578125}, {0.3515625, -0.0234375, 0.6171875},
        {-0.3515625, -0.0234375, 0.6171875}, {0.3515625, 0.03125, 0.71875},
        {-0.3515625, 0.03125, 0.71875}, {0.3515625, 0.1328125, 0.78125},
        {-0.3515625, 0.1328125, 0.78125}, {0.2734375, 0.1640625, 0.796875},
        {-0.2734375, 0.1640625, 0.796875}, {0.203125, 0.09375, 0.7421875},
        {-0.203125, 0.09375, 0.7421875}, {0.15625, 0.0546875, 0.6484375},
        {-0.15625, 0.0546875, 0.6484375}, {0.078125, 0.2421875, 0.65625},
        {-0.078125, 0.2421875, 0.65625}, {0.140625, 0.2421875, 0.7421875},
        {-0.140625, 0.2421875, 0.7421875}, {0.2421875, 0.2421875, 0.796875},
        {-0.2421875, 0.2421875, 0.796875}, {0.2734375, 0.328125, 0.796875},
        {-0.2734375, 0.328125, 0.796875}, {0.203125, 0.390625, 0.7421875},
        {-0.203125, 0.390625, 0.7421875}, {0.15625, 0.4375, 0.6484375},
        {-0.15625, 0.4375, 0.6484375}, {0.3515625, 0.515625, 0.6171875},
        {-0.3515625, 0.515625, 0.6171875}, {0.3515625, 0.453125, 0.71875},
        {-0.3515625, 0.453125, 0.71875}, {0.3515625, 0.359375, 0.78125},
        {-0.3515625, 0.359375, 0.78125}, {0.4375, 0.328125, 0.765625},
        {-0.4375, 0.328125, 0.765625}, {0.5, 0.390625, 0.6875},
        {-0.5, 0.390625, 0.6875}, {0.546875, 0.4375, 0.578125},
        {-0.546875, 0.4375, 0.578125}, {0.625, 0.2421875, 0.5625},
        {-0.625, 0.2421875, 0.5625}, {0.5625, 0.2421875, 0.671875},
        {-0.5625, 0.2421875, 0.671875}, {0.46875, 0.2421875, 0.7578125},
        {-0.46875, 0.2421875, 0.7578125}, {0.4765625, 0.2421875, 0.7734375},
        {-0.4765625, 0.2421875, 0.7734375}, {0.4453125, 0.3359375, 0.78125},
        {-0.4453125, 0.3359375, 0.78125}, {0.3515625, 0.375, 0.8046875},
        {-0.3515625, 0.375, 0.8046875}, {0.265625, 0.3359375, 0.8203125},
        {-0.265625, 0.3359375, 0.8203125}, {0.2265625, 0.2421875, 0.8203125},
        {-0.2265625, 0.2421875, 0.8203125}, {0.265625, 0.15625, 0.8203125},
        {-0.265625, 0.15625, 0.8203125}, {0.3515625, 0.2421875, 0.828125},
        {-0.3515625, 0.2421875, 0.828125}, {0.3515625, 0.1171875, 0.8046875},
        {-0.3515625, 0.1171875, 0.8046875}, {0.4453125, 0.15625, 0.78125},
        {-0.4453125, 0.15625, 0.78125}, {0.0, 0.4296875, 0.7421875},
        {0.0, 0.3515625, 0.8203125}, {0.0, -0.6796875, 0.734375},
        {0.0, -0.3203125, 0.78125}, {0.0, -0.1875, 0.796875},
        {0.0, -0.7734375, 0.71875}, {0.0, 0.40625, 0.6015625},
        {0.0, 0.5703125, 0.5703125}, {0.0, 0.8984375, -0.546875},
        {0.0, 0.5625, -0.8515625}, {0.0, 0.0703125, -0.828125},
        {0.0, -0.3828125, -0.3515625}, {0.203125, -0.1875, 0.5625},
        {-0.203125, -0.1875, 0.5625}, {0.3125, -0.4375, 0.5703125},
        {-0.3125, -0.4375, 0.5703125}, {0.3515625, -0.6953125, 0.5703125},
        {-0.3515625, -0.6953125, 0.5703125}, {0.3671875, -0.890625, 0.53125},
        {-0.3671875, -0.890625, 0.53125}, {0.328125, -0.9453125, 0.5234375},
        {-0.328125, -0.9453125, 0.5234375}, {0.1796875, -0.96875, 0.5546875},
        {-0.1796875, -0.96875, 0.5546875}, {0.0, -0.984375, 0.578125},
        {0.4375, -0.140625, 0.53125}, {-0.4375, -0.140625, 0.53125},
        {0.6328125, -0.0390625, 0.5390625}, {-0.6328125, -0.0390625, 0.5390625},
        {0.828125, 0.1484375, 0.4453125}, {-0.828125, 0.1484375, 0.4453125},
        {0.859375, 0.4296875, 0.59375}, {-0.859375, 0.4296875, 0.59375},
        {0.7109375, 0.484375, 0.625}, {-0.7109375, 0.484375, 0.625},
        {0.4921875, 0.6015625, 0.6875}, {-0.4921875, 0.6015625, 0.6875},
        {0.3203125, 0.7578125, 0.734375}, {-0.3203125, 0.7578125, 0.734375},
        {0.15625, 0.71875, 0.7578125}, {-0.15625, 0.71875, 0.7578125},
        {0.0625, 0.4921875, 0.75}, {-0.0625, 0.4921875, 0.75},
        {0.1640625, 0.4140625, 0.7734375}, {-0.1640625, 0.4140625, 0.7734375},
        {0.125, 0.3046875, 0.765625}, {-0.125, 0.3046875, 0.765625},
        {0.203125, 0.09375, 0.7421875}, {-0.203125, 0.09375, 0.7421875},
        {0.375, 0.015625, 0.703125}, {-0.375, 0.015625, 0.703125},
        {0.4921875, 0.0625, 0.671875}, {-0.4921875, 0.0625, 0.671875},
        {0.625, 0.1875, 0.6484375}, {-0.625, 0.1875, 0.6484375},
        {0.640625, 0.296875, 0.6484375}, {-0.640625, 0.296875, 0.6484375},
        {0.6015625, 0.375, 0.6640625}, {-0.6015625, 0.375, 0.6640625},
        {0.4296875, 0.4375, 0.71875}, {-0.4296875, 0.4375, 0.71875},
        {0.25, 0.46875, 0.7578125}, {-0.25, 0.46875, 0.7578125},
        {0.0, -0.765625, 0.734375}, {0.109375, -0.71875, 0.734375},
        {-0.109375, -0.71875, 0.734375}, {0.1171875, -0.8359375, 0.7109375},
        {-0.1171875, -0.8359375, 0.7109375}, {0.0625, -0.8828125, 0.6953125},
        {-0.0625, -0.8828125, 0.6953125}, {0.0, -0.890625, 0.6875},
        {0.0, -0.1953125, 0.75}, {0.0, -0.140625, 0.7421875},
        {0.1015625, -0.1484375, 0.7421875}, {-0.1015625, -0.1484375, 0.7421875},
        {0.125, -0.2265625, 0.75}, {-0.125, -0.2265625, 0.75},
        {0.0859375, -0.2890625, 0.7421875}, {-0.0859375, -0.2890625, 0.7421875},
        {0.3984375, -0.046875, 0.671875}, {-0.3984375, -0.046875, 0.671875},
        {0.6171875, 0.0546875, 0.625}, {-0.6171875, 0.0546875, 0.625},
        {0.7265625, 0.203125, 0.6015625}, {-0.7265625, 0.203125, 0.6015625},
        {0.7421875, 0.375, 0.65625}, {-0.7421875, 0.375, 0.65625},
        {0.6875, 0.4140625, 0.7265625}, {-0.6875, 0.4140625, 0.7265625},
        {0.4375, 0.546875, 0.796875}, {-0.4375, 0.546875, 0.796875},
        {0.3125, 0.640625, 0.8359375}, {-0.3125, 0.640625, 0.8359375},
        {0.203125, 0.6171875, 0.8515625}, {-0.203125, 0.6171875, 0.8515625},
        {0.1015625, 0.4296875, 0.84375}, {-0.1015625, 0.4296875, 0.84375},
        {0.125, -0.1015625, 0.8125}, {-0.125, -0.1015625, 0.8125},
        {0.2109375, -0.4453125, 0.7109375}, {-0.2109375, -0.4453125, 0.7109375},
        {0.25, -0.703125, 0.6875}, {-0.25, -0.703125, 0.6875},
        {0.265625, -0.8203125, 0.6640625}, {-0.265625, -0.8203125, 0.6640625},
        {0.234375, -0.9140625, 0.6328125}, {-0.234375, -0.9140625, 0.6328125},
        {0.1640625, -0.9296875, 0.6328125}, {-0.1640625, -0.9296875, 0.6328125},
        {0.0, -0.9453125, 0.640625}, {0.0, 0.046875, 0.7265625},
        {0.0, 0.2109375, 0.765625}, {0.328125, 0.4765625, 0.7421875},
        {-0.328125, 0.4765625, 0.7421875}, {0.1640625, 0.140625, 0.75},
        {-0.1640625, 0.140625, 0.75}, {0.1328125, 0.2109375, 0.7578125},
        {-0.1328125, 0.2109375, 0.7578125}, {0.1171875, -0.6875, 0.734375},
        {-0.1171875, -0.6875, 0.734375}, {0.078125, -0.4453125, 0.75},
        {-0.078125, -0.4453125, 0.75}, {0.0, -0.4453125, 0.75},
        {0.0, -0.328125, 0.7421875}, {0.09375, -0.2734375, 0.78125},
        {-0.09375, -0.2734375, 0.78125}, {0.1328125, -0.2265625, 0.796875},
        {-0.1328125, -0.2265625, 0.796875}, {0.109375, -0.1328125, 0.78125},
        {-0.109375, -0.1328125, 0.78125}, {0.0390625, -0.125, 0.78125},
        {-0.0390625, -0.125, 0.78125}, {0.0, -0.203125, 0.828125},
        {0.046875, -0.1484375, 0.8125}, {-0.046875, -0.1484375, 0.8125},
        {0.09375, -0.15625, 0.8125}, {-0.09375, -0.15625, 0.8125},
        {0.109375, -0.2265625, 0.828125}, {-0.109375, -0.2265625, 0.828125},
        {0.078125, -0.25, 0.8046875}, {-0.078125, -0.25, 0.8046875},
        {0.0, -0.2890625, 0.8046875}, {0.2578125, -0.3125, 0.5546875},
        {-0.2578125, -0.3125, 0.5546875}, {0.1640625, -0.2421875, 0.7109375},
        {-0.1640625, -0.2421875, 0.7109375}, {0.1796875, -0.3125, 0.7109375},
        {-0.1796875, -0.3125, 0.7109375}, {0.234375, -0.25, 0.5546875},
        {-0.234375, -0.25, 0.5546875}, {0.0, -0.875, 0.6875},
        {0.046875, -0.8671875, 0.6875}, {-0.046875, -0.8671875, 0.6875},
        {0.09375, -0.8203125, 0.7109375}, {-0.09375, -0.8203125, 0.7109375},
        {0.09375, -0.7421875, 0.7265625}, {-0.09375, -0.7421875, 0.7265625},
        {0.0, -0.78125, 0.65625}, {0.09375, -0.75, 0.6640625},
        {-0.09375, -0.75, 0.6640625}, {0.09375, -0.8125, 0.640625},
        {-0.09375, -0.8125, 0.640625}, {0.046875, -0.8515625, 0.6328125},
        {-0.046875, -0.8515625, 0.6328125}, {0.0, -0.859375, 0.6328125},
        {0.171875, 0.21875, 0.78125}, {-0.171875, 0.21875, 0.78125},
        {0.1875, 0.15625, 0.7734375}, {-0.1875, 0.15625, 0.7734375},
        {0.3359375, 0.4296875, 0.7578125}, {-0.3359375, 0.4296875, 0.7578125},
        {0.2734375, 0.421875, 0.7734375}, {-0.2734375, 0.421875, 0.7734375},
        {0.421875, 0.3984375, 0.7734375}, {-0.421875, 0.3984375, 0.7734375},
        {0.5625, 0.3515625, 0.6953125}, {-0.5625, 0.3515625, 0.6953125},
        {0.5859375, 0.2890625, 0.6875}, {-0.5859375, 0.2890625, 0.6875},
        {0.578125, 0.1953125, 0.6796875}, {-0.578125, 0.1953125, 0.6796875},
        {0.4765625, 0.1015625, 0.71875}, {-0.4765625, 0.1015625, 0.71875},
        {0.375, 0.0625, 0.7421875}, {-0.375, 0.0625, 0.7421875},
        {0.2265625, 0.109375, 0.78125}, {-0.2265625, 0.109375, 0.78125},
        {0.1796875, 0.296875, 0.78125}, {-0.1796875, 0.296875, 0.78125},
        {0.2109375, 0.375, 0.78125}, {-0.2109375, 0.375, 0.78125},
        {0.234375, 0.359375, 0.7578125}, {-0.234375, 0.359375, 0.7578125},
        {0.1953125, 0.296875, 0.7578125}, {-0.1953125, 0.296875, 0.7578125},
        {0.2421875, 0.125, 0.7578125}, {-0.2421875, 0.125, 0.7578125},
        {0.375, 0.0859375, 0.7265625}, {-0.375, 0.0859375, 0.7265625},
        {0.4609375, 0.1171875, 0.703125}, {-0.4609375, 0.1171875, 0.703125},
        {0.546875, 0.2109375, 0.671875}, {-0.546875, 0.2109375, 0.671875},
        {0.5546875, 0.28125, 0.671875}, {-0.5546875, 0.28125, 0.671875},
        {0.53125, 0.3359375, 0.6796875}, {-0.53125, 0.3359375, 0.6796875},
        {0.4140625, 0.390625, 0.75}, {-0.4140625, 0.390625, 0.75},
        {0.28125, 0.3984375, 0.765625}, {-0.28125, 0.3984375, 0.765625},
        {0.3359375, 0.40625, 0.75}, {-0.3359375, 0.40625, 0.75},
        {0.203125, 0.171875, 0.75}, {-0.203125, 0.171875, 0.75},
        {0.1953125, 0.2265625, 0.75}, {-0.1953125, 0.2265625, 0.75},
        {0.109375, 0.4609375, 0.609375}, {-0.109375, 0.4609375, 0.609375},
        {0.1953125, 0.6640625, 0.6171875}, {-0.1953125, 0.6640625, 0.6171875},
        {0.3359375, 0.6875, 0.59375}, {-0.3359375, 0.6875, 0.59375},
        {0.484375, 0.5546875, 0.5546875}, {-0.484375, 0.5546875, 0.5546875},
        {0.6796875, 0.453125, 0.4921875}, {-0.6796875, 0.453125, 0.4921875},
        {0.796875, 0.40625, 0.4609375}, {-0.796875, 0.40625, 0.4609375},
        {0.7734375, 0.1640625, 0.375}, {-0.7734375, 0.1640625, 0.375},
        {0.6015625, 0.0, 0.4140625}, {-0.6015625, 0.0, 0.4140625},
        {0.4375, -0.09375, 0.46875}, {-0.4375, -0.09375, 0.46875},
        {0.0, 0.8984375, 0.2890625}, {0.0, 0.984375, -0.078125},
        {0.0, -0.1953125, -0.671875}, {0.0, -0.4609375, 0.1875},
        {0.0, -0.9765625, 0.4609375}, {0.0, -0.8046875, 0.34375},
        {0.0, -0.5703125, 0.3203125}, {0.0, -0.484375, 0.28125},
        {0.8515625, 0.234375, 0.0546875}, {-0.8515625, 0.234375, 0.0546875},
        {0.859375, 0.3203125, -0.046875}, {-0.859375, 0.3203125, -0.046875},
        {0.7734375, 0.265625, -0.4375}, {-0.7734375, 0.265625, -0.4375},
        {0.4609375, 0.4375, -0.703125}, {-0.4609375, 0.4375, -0.703125},
        {0.734375, -0.046875, 0.0703125}, {-0.734375, -0.046875, 0.0703125},
        {0.59375, -0.125, -0.1640625}, {-0.59375, -0.125, -0.1640625},
        {0.640625, -0.0078125, -0.4296875}, {-0.640625, -0.0078125, -0.4296875},
        {0.3359375, 0.0546875, -0.6640625}, {-0.3359375, 0.0546875, -0.6640625},
        {0.234375, -0.3515625, 0.40625}, {-0.234375, -0.3515625, 0.40625},
        {0.1796875, -0.4140625, 0.2578125}, {-0.1796875, -0.4140625, 0.2578125},
        {0.2890625, -0.7109375, 0.3828125}, {-0.2890625, -0.7109375, 0.3828125},
        {0.25, -0.5, 0.390625}, {-0.25, -0.5, 0.390625},
        {0.328125, -0.9140625, 0.3984375}, {-0.328125, -0.9140625, 0.3984375},
        {0.140625, -0.7578125, 0.3671875}, {-0.140625, -0.7578125, 0.3671875},
        {0.125, -0.5390625, 0.359375}, {-0.125, -0.5390625, 0.359375},
        {0.1640625, -0.9453125, 0.4375}, {-0.1640625, -0.9453125, 0.4375},
        {0.21875, -0.28125, 0.4296875}, {-0.21875, -0.28125, 0.4296875},
        {0.2109375, -0.2265625, 0.46875}, {-0.2109375, -0.2265625, 0.46875},
        {0.203125, -0.171875, 0.5}, {-0.203125, -0.171875, 0.5},
        {0.2109375, -0.390625, 0.1640625}, {-0.2109375, -0.390625, 0.1640625},
        {0.296875, -0.3125, -0.265625}, {-0.296875, -0.3125, -0.265625},
        {0.34375, -0.1484375, -0.5390625}, {-0.34375, -0.1484375, -0.5390625},
        {0.453125, 0.8671875, -0.3828125}, {-0.453125, 0.8671875, -0.3828125},
        {0.453125, 0.9296875, -0.0703125}, {-0.453125, 0.9296875, -0.0703125},
        {0.453125, 0.8515625, 0.234375}, {-0.453125, 0.8515625, 0.234375},
        {0.4609375, 0.5234375, 0.4296875}, {-0.4609375, 0.5234375, 0.4296875},
        {0.7265625, 0.40625, 0.3359375}, {-0.7265625, 0.40625, 0.3359375},
        {0.6328125, 0.453125, 0.28125}, {-0.6328125, 0.453125, 0.28125},
        {0.640625, 0.703125, 0.0546875}, {-0.640625, 0.703125, 0.0546875},
        {0.796875, 0.5625, 0.125}, {-0.796875, 0.5625, 0.125},
        {0.796875, 0.6171875, -0.1171875}, {-0.796875, 0.6171875, -0.1171875},
        {0.640625, 0.75, -0.1953125}, {-0.640625, 0.75, -0.1953125},
        {0.640625, 0.6796875, -0.4453125}, {-0.640625, 0.6796875, -0.4453125},
        {0.796875, 0.5390625, -0.359375}, {-0.796875, 0.5390625, -0.359375},
        {0.6171875, 0.328125, -0.5859375}, {-0.6171875, 0.328125, -0.5859375},
        {0.484375, 0.0234375, -0.546875}, {-0.484375, 0.0234375, -0.546875},
        {0.8203125, 0.328125, -0.203125}, {-0.8203125, 0.328125, -0.203125},
        {0.40625, -0.171875, 0.1484375}, {-0.40625, -0.171875, 0.1484375},
        {0.4296875, -0.1953125, -0.2109375},
        {-0.4296875, -0.1953125, -0.2109375}, {0.890625, 0.40625, -0.234375},
        {-0.890625, 0.40625, -0.234375}, {0.7734375, -0.140625, -0.125},
        {-0.7734375, -0.140625, -0.125}, {1.0390625, -0.1015625, -0.328125},
        {-1.0390625, -0.1015625, -0.328125}, {1.28125, 0.0546875, -0.4296875},
        {-1.28125, 0.0546875, -0.4296875}, {1.3515625, 0.3203125, -0.421875},
        {-1.3515625, 0.3203125, -0.421875}, {1.234375, 0.5078125, -0.421875},
        {-1.234375, 0.5078125, -0.421875}, {1.0234375, 0.4765625, -0.3125},
        {-1.0234375, 0.4765625, -0.3125}, {1.015625, 0.4140625, -0.2890625},
        {-1.015625, 0.4140625, -0.2890625}, {1.1875, 0.4375, -0.390625},
        {-1.1875, 0.4375, -0.390625}, {1.265625, 0.2890625, -0.40625},
        {-1.265625, 0.2890625, -0.40625}, {1.2109375, 0.078125, -0.40625},
        {-1.2109375, 0.078125, -0.40625}, {1.03125, -0.0390625, -0.3046875},
        {-1.03125, -0.0390625, -0.3046875}, {0.828125, -0.0703125, -0.1328125},
        {-0.828125, -0.0703125, -0.1328125}, {0.921875, 0.359375, -0.21875},
        {-0.921875, 0.359375, -0.21875}, {0.9453125, 0.3046875, -0.2890625},
        {-0.9453125, 0.3046875, -0.2890625},
        {0.8828125, -0.0234375, -0.2109375},
        {-0.8828125, -0.0234375, -0.2109375}, {1.0390625, 0.0, -0.3671875},
        {-1.0390625, 0.0, -0.3671875}, {1.1875, 0.09375, -0.4453125},
        {-1.1875, 0.09375, -0.4453125}, {1.234375, 0.25, -0.4453125},
        {-1.234375, 0.25, -0.4453125}, {1.171875, 0.359375, -0.4375},
        {-1.171875, 0.359375, -0.4375}, {1.0234375, 0.34375, -0.359375},
        {-1.0234375, 0.34375, -0.359375}, {0.84375, 0.2890625, -0.2109375},
        {-0.84375, 0.2890625, -0.2109375}, {0.8359375, 0.171875, -0.2734375},
        {-0.8359375, 0.171875, -0.2734375}, {0.7578125, 0.09375, -0.2734375},
        {-0.7578125, 0.09375, -0.2734375}, {0.8203125, 0.0859375, -0.2734375},
        {-0.8203125, 0.0859375, -0.2734375}, {0.84375, 0.015625, -0.2734375},
        {-0.84375, 0.015625, -0.2734375}, {0.8125, -0.015625, -0.2734375},
        {-0.8125, -0.015625, -0.2734375}, {0.7265625, 0.0, -0.0703125},
        {-0.7265625, 0.0, -0.0703125}, {0.71875, -0.0234375, -0.171875},
        {-0.71875, -0.0234375, -0.171875}, {0.71875, 0.0390625, -0.1875},
        {-0.71875, 0.0390625, -0.1875}, {0.796875, 0.203125, -0.2109375},
        {-0.796875, 0.203125, -0.2109375}, {0.890625, 0.2421875, -0.265625},
        {-0.890625, 0.2421875, -0.265625}, {0.890625, 0.234375, -0.3203125},
        {-0.890625, 0.234375, -0.3203125}, {0.8125, -0.015625, -0.3203125},
        {-0.8125, -0.015625, -0.3203125}, {0.8515625, 0.015625, -0.3203125},
        {-0.8515625, 0.015625, -0.3203125}, {0.828125, 0.078125, -0.3203125},
        {-0.828125, 0.078125, -0.3203125}, {0.765625, 0.09375, -0.3203125},
        {-0.765625, 0.09375, -0.3203125}, {0.84375, 0.171875, -0.3203125},
        {-0.84375, 0.171875, -0.3203125}, {1.0390625, 0.328125, -0.4140625},
        {-1.0390625, 0.328125, -0.4140625}, {1.1875, 0.34375, -0.484375},
        {-1.1875, 0.34375, -0.484375}, {1.2578125, 0.2421875, -0.4921875},
        {-1.2578125, 0.2421875, -0.4921875}, {1.2109375, 0.0859375, -0.484375},
        {-1.2109375, 0.0859375, -0.484375}, {1.046875, 0.0, -0.421875},
        {-1.046875, 0.0, -0.421875}, {0.8828125, -0.015625, -0.265625},
        {-0.8828125, -0.015625, -0.265625}, {0.953125, 0.2890625, -0.34375},
        {-0.953125, 0.2890625, -0.34375}, {0.890625, 0.109375, -0.328125},
        {-0.890625, 0.109375, -0.328125}, {0.9375, 0.0625, -0.3359375},
        {-0.9375, 0.0625, -0.3359375}, {1.0, 0.125, -0.3671875},
        {-1.0, 0.125, -0.3671875}, {0.9609375, 0.171875, -0.3515625},
        {-0.9609375, 0.171875, -0.3515625}, {1.015625, 0.234375, -0.375},
        {-1.015625, 0.234375, -0.375}, {1.0546875, 0.1875, -0.3828125},
        {-1.0546875, 0.1875, -0.3828125}, {1.109375, 0.2109375, -0.390625},
        {-1.109375, 0.2109375, -0.390625}, {1.0859375, 0.2734375, -0.390625},
        {-1.0859375, 0.2734375, -0.390625}, {1.0234375, 0.4375, -0.484375},
        {-1.0234375, 0.4375, -0.484375}, {1.25, 0.46875, -0.546875},
        {-1.25, 0.46875, -0.546875}, {1.3671875, 0.296875, -0.5},
        {-1.3671875, 0.296875, -0.5}, {1.3125, 0.0546875, -0.53125},
        {-1.3125, 0.0546875, -0.53125}, {1.0390625, -0.0859375, -0.4921875},
        {-1.0390625, -0.0859375, -0.4921875}, {0.7890625, -0.125, -0.328125},
        {-0.7890625, -0.125, -0.328125}, {0.859375, 0.3828125, -0.3828125},
        {-0.859375, 0.3828125, -0.3828125}};
    static auto suzanne_triangles = std::vector<vec3i>{{60, 64, 48},
        {49, 65, 61}, {62, 64, 60}, {61, 65, 63}, {60, 58, 62}, {63, 59, 61},
        {60, 56, 58}, {59, 57, 61}, {60, 54, 56}, {57, 55, 61}, {60, 52, 54},
        {55, 53, 61}, {60, 50, 52}, {53, 51, 61}, {60, 48, 50}, {51, 49, 61},
        {224, 228, 226}, {227, 229, 225}, {72, 283, 73}, {73, 284, 72},
        {341, 347, 383}, {384, 348, 342}, {299, 345, 343}, {344, 346, 300},
        {323, 379, 351}, {352, 380, 324}, {441, 443, 445}, {446, 444, 442},
        {463, 491, 465}, {466, 492, 464}, {495, 497, 499}, {500, 498, 496}};
    static auto suzanne_quads = std::vector<vec4i>{{46, 0, 2, 44},
        {3, 1, 47, 45}, {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4},
        {7, 9, 3, 5}, {0, 10, 8, 2}, {9, 11, 1, 3}, {10, 12, 14, 8},
        {15, 13, 11, 9}, {8, 14, 16, 6}, {17, 15, 9, 7}, {14, 20, 18, 16},
        {19, 21, 15, 17}, {12, 22, 20, 14}, {21, 23, 13, 15}, {22, 24, 26, 20},
        {27, 25, 23, 21}, {20, 26, 28, 18}, {29, 27, 21, 19}, {26, 32, 30, 28},
        {31, 33, 27, 29}, {24, 34, 32, 26}, {33, 35, 25, 27}, {34, 36, 38, 32},
        {39, 37, 35, 33}, {32, 38, 40, 30}, {41, 39, 33, 31}, {38, 44, 42, 40},
        {43, 45, 39, 41}, {36, 46, 44, 38}, {45, 47, 37, 39}, {46, 36, 50, 48},
        {51, 37, 47, 49}, {36, 34, 52, 50}, {53, 35, 37, 51}, {34, 24, 54, 52},
        {55, 25, 35, 53}, {24, 22, 56, 54}, {57, 23, 25, 55}, {22, 12, 58, 56},
        {59, 13, 23, 57}, {12, 10, 62, 58}, {63, 11, 13, 59}, {10, 0, 64, 62},
        {65, 1, 11, 63}, {0, 46, 48, 64}, {49, 47, 1, 65}, {88, 173, 175, 90},
        {175, 174, 89, 90}, {86, 171, 173, 88}, {174, 172, 87, 89},
        {84, 169, 171, 86}, {172, 170, 85, 87}, {82, 167, 169, 84},
        {170, 168, 83, 85}, {80, 165, 167, 82}, {168, 166, 81, 83},
        {78, 91, 145, 163}, {146, 92, 79, 164}, {91, 93, 147, 145},
        {148, 94, 92, 146}, {93, 95, 149, 147}, {150, 96, 94, 148},
        {95, 97, 151, 149}, {152, 98, 96, 150}, {97, 99, 153, 151},
        {154, 100, 98, 152}, {99, 101, 155, 153}, {156, 102, 100, 154},
        {101, 103, 157, 155}, {158, 104, 102, 156}, {103, 105, 159, 157},
        {160, 106, 104, 158}, {105, 107, 161, 159}, {162, 108, 106, 160},
        {107, 66, 67, 161}, {67, 66, 108, 162}, {109, 127, 159, 161},
        {160, 128, 110, 162}, {127, 178, 157, 159}, {158, 179, 128, 160},
        {125, 155, 157, 178}, {158, 156, 126, 179}, {123, 153, 155, 125},
        {156, 154, 124, 126}, {121, 151, 153, 123}, {154, 152, 122, 124},
        {119, 149, 151, 121}, {152, 150, 120, 122}, {117, 147, 149, 119},
        {150, 148, 118, 120}, {115, 145, 147, 117}, {148, 146, 116, 118},
        {113, 163, 145, 115}, {146, 164, 114, 116}, {113, 180, 176, 163},
        {176, 181, 114, 164}, {109, 161, 67, 111}, {67, 162, 110, 112},
        {111, 67, 177, 182}, {177, 67, 112, 183}, {176, 180, 182, 177},
        {183, 181, 176, 177}, {134, 136, 175, 173}, {175, 136, 135, 174},
        {132, 134, 173, 171}, {174, 135, 133, 172}, {130, 132, 171, 169},
        {172, 133, 131, 170}, {165, 186, 184, 167}, {185, 187, 166, 168},
        {130, 169, 167, 184}, {168, 170, 131, 185}, {143, 189, 188, 186},
        {188, 189, 144, 187}, {184, 186, 188, 68}, {188, 187, 185, 68},
        {129, 130, 184, 68}, {185, 131, 129, 68}, {141, 192, 190, 143},
        {191, 193, 142, 144}, {139, 194, 192, 141}, {193, 195, 140, 142},
        {138, 196, 194, 139}, {195, 197, 138, 140}, {137, 70, 196, 138},
        {197, 70, 137, 138}, {189, 143, 190, 69}, {191, 144, 189, 69},
        {69, 190, 205, 207}, {206, 191, 69, 207}, {70, 198, 199, 196},
        {200, 198, 70, 197}, {196, 199, 201, 194}, {202, 200, 197, 195},
        {194, 201, 203, 192}, {204, 202, 195, 193}, {192, 203, 205, 190},
        {206, 204, 193, 191}, {198, 203, 201, 199}, {202, 204, 198, 200},
        {198, 207, 205, 203}, {206, 207, 198, 204}, {138, 139, 163, 176},
        {164, 140, 138, 176}, {139, 141, 210, 163}, {211, 142, 140, 164},
        {141, 143, 212, 210}, {213, 144, 142, 211}, {143, 186, 165, 212},
        {166, 187, 144, 213}, {80, 208, 212, 165}, {213, 209, 81, 166},
        {208, 214, 210, 212}, {211, 215, 209, 213}, {78, 163, 210, 214},
        {211, 164, 79, 215}, {130, 129, 71, 221}, {71, 129, 131, 222},
        {132, 130, 221, 219}, {222, 131, 133, 220}, {134, 132, 219, 217},
        {220, 133, 135, 218}, {136, 134, 217, 216}, {218, 135, 136, 216},
        {216, 217, 228, 230}, {229, 218, 216, 230}, {217, 219, 226, 228},
        {227, 220, 218, 229}, {219, 221, 224, 226}, {225, 222, 220, 227},
        {221, 71, 223, 224}, {223, 71, 222, 225}, {223, 230, 228, 224},
        {229, 230, 223, 225}, {182, 180, 233, 231}, {234, 181, 183, 232},
        {111, 182, 231, 253}, {232, 183, 112, 254}, {109, 111, 253, 255},
        {254, 112, 110, 256}, {180, 113, 251, 233}, {252, 114, 181, 234},
        {113, 115, 249, 251}, {250, 116, 114, 252}, {115, 117, 247, 249},
        {248, 118, 116, 250}, {117, 119, 245, 247}, {246, 120, 118, 248},
        {119, 121, 243, 245}, {244, 122, 120, 246}, {121, 123, 241, 243},
        {242, 124, 122, 244}, {123, 125, 239, 241}, {240, 126, 124, 242},
        {125, 178, 235, 239}, {236, 179, 126, 240}, {178, 127, 237, 235},
        {238, 128, 179, 236}, {127, 109, 255, 237}, {256, 110, 128, 238},
        {237, 255, 257, 275}, {258, 256, 238, 276}, {235, 237, 275, 277},
        {276, 238, 236, 278}, {239, 235, 277, 273}, {278, 236, 240, 274},
        {241, 239, 273, 271}, {274, 240, 242, 272}, {243, 241, 271, 269},
        {272, 242, 244, 270}, {245, 243, 269, 267}, {270, 244, 246, 268},
        {247, 245, 267, 265}, {268, 246, 248, 266}, {249, 247, 265, 263},
        {266, 248, 250, 264}, {251, 249, 263, 261}, {264, 250, 252, 262},
        {233, 251, 261, 279}, {262, 252, 234, 280}, {255, 253, 259, 257},
        {260, 254, 256, 258}, {253, 231, 281, 259}, {282, 232, 254, 260},
        {231, 233, 279, 281}, {280, 234, 232, 282}, {66, 107, 283, 72},
        {284, 108, 66, 72}, {107, 105, 285, 283}, {286, 106, 108, 284},
        {105, 103, 287, 285}, {288, 104, 106, 286}, {103, 101, 289, 287},
        {290, 102, 104, 288}, {101, 99, 291, 289}, {292, 100, 102, 290},
        {99, 97, 293, 291}, {294, 98, 100, 292}, {97, 95, 295, 293},
        {296, 96, 98, 294}, {95, 93, 297, 295}, {298, 94, 96, 296},
        {93, 91, 299, 297}, {300, 92, 94, 298}, {307, 308, 327, 337},
        {328, 308, 307, 338}, {306, 307, 337, 335}, {338, 307, 306, 336},
        {305, 306, 335, 339}, {336, 306, 305, 340}, {88, 90, 305, 339},
        {305, 90, 89, 340}, {86, 88, 339, 333}, {340, 89, 87, 334},
        {84, 86, 333, 329}, {334, 87, 85, 330}, {82, 84, 329, 331},
        {330, 85, 83, 332}, {329, 335, 337, 331}, {338, 336, 330, 332},
        {329, 333, 339, 335}, {340, 334, 330, 336}, {325, 331, 337, 327},
        {338, 332, 326, 328}, {80, 82, 331, 325}, {332, 83, 81, 326},
        {208, 341, 343, 214}, {344, 342, 209, 215}, {80, 325, 341, 208},
        {342, 326, 81, 209}, {78, 214, 343, 345}, {344, 215, 79, 346},
        {78, 345, 299, 91}, {300, 346, 79, 92}, {76, 323, 351, 303},
        {352, 324, 76, 303}, {303, 351, 349, 77}, {350, 352, 303, 77},
        {77, 349, 347, 304}, {348, 350, 77, 304}, {304, 347, 327, 308},
        {328, 348, 304, 308}, {325, 327, 347, 341}, {348, 328, 326, 342},
        {295, 297, 317, 309}, {318, 298, 296, 310}, {75, 315, 323, 76},
        {324, 316, 75, 76}, {301, 357, 355, 302}, {356, 358, 301, 302},
        {302, 355, 353, 74}, {354, 356, 302, 74}, {74, 353, 315, 75},
        {316, 354, 74, 75}, {291, 293, 361, 363}, {362, 294, 292, 364},
        {363, 361, 367, 365}, {368, 362, 364, 366}, {365, 367, 369, 371},
        {370, 368, 366, 372}, {371, 369, 375, 373}, {376, 370, 372, 374},
        {313, 377, 373, 375}, {374, 378, 314, 376}, {315, 353, 373, 377},
        {374, 354, 316, 378}, {353, 355, 371, 373}, {372, 356, 354, 374},
        {355, 357, 365, 371}, {366, 358, 356, 372}, {357, 359, 363, 365},
        {364, 360, 358, 366}, {289, 291, 363, 359}, {364, 292, 290, 360},
        {73, 359, 357, 301}, {358, 360, 73, 301}, {283, 285, 287, 289},
        {288, 286, 284, 290}, {283, 289, 359, 73}, {360, 290, 284, 73},
        {293, 295, 309, 361}, {310, 296, 294, 362}, {309, 311, 367, 361},
        {368, 312, 310, 362}, {311, 381, 369, 367}, {370, 382, 312, 368},
        {313, 375, 369, 381}, {370, 376, 314, 382}, {347, 349, 385, 383},
        {386, 350, 348, 384}, {317, 383, 385, 319}, {386, 384, 318, 320},
        {297, 299, 383, 317}, {384, 300, 298, 318}, {299, 343, 341, 383},
        {342, 344, 300, 384}, {313, 321, 379, 377}, {380, 322, 314, 378},
        {315, 377, 379, 323}, {380, 378, 316, 324}, {319, 385, 379, 321},
        {380, 386, 320, 322}, {349, 351, 379, 385}, {380, 352, 350, 386},
        {399, 387, 413, 401}, {414, 388, 400, 402}, {399, 401, 403, 397},
        {404, 402, 400, 398}, {397, 403, 405, 395}, {406, 404, 398, 396},
        {395, 405, 407, 393}, {408, 406, 396, 394}, {393, 407, 409, 391},
        {410, 408, 394, 392}, {391, 409, 411, 389}, {412, 410, 392, 390},
        {409, 419, 417, 411}, {418, 420, 410, 412}, {407, 421, 419, 409},
        {420, 422, 408, 410}, {405, 423, 421, 407}, {422, 424, 406, 408},
        {403, 425, 423, 405}, {424, 426, 404, 406}, {401, 427, 425, 403},
        {426, 428, 402, 404}, {401, 413, 415, 427}, {416, 414, 402, 428},
        {317, 319, 443, 441}, {444, 320, 318, 442}, {319, 389, 411, 443},
        {412, 390, 320, 444}, {309, 317, 441, 311}, {442, 318, 310, 312},
        {381, 429, 413, 387}, {414, 430, 382, 388}, {411, 417, 439, 443},
        {440, 418, 412, 444}, {437, 445, 443, 439}, {444, 446, 438, 440},
        {433, 445, 437, 435}, {438, 446, 434, 436}, {431, 447, 445, 433},
        {446, 448, 432, 434}, {429, 447, 431, 449}, {432, 448, 430, 450},
        {413, 429, 449, 415}, {450, 430, 414, 416}, {311, 447, 429, 381},
        {430, 448, 312, 382}, {311, 441, 445, 447}, {446, 442, 312, 448},
        {415, 449, 451, 475}, {452, 450, 416, 476}, {449, 431, 461, 451},
        {462, 432, 450, 452}, {431, 433, 459, 461}, {460, 434, 432, 462},
        {433, 435, 457, 459}, {458, 436, 434, 460}, {435, 437, 455, 457},
        {456, 438, 436, 458}, {437, 439, 453, 455}, {454, 440, 438, 456},
        {439, 417, 473, 453}, {474, 418, 440, 454}, {427, 415, 475, 463},
        {476, 416, 428, 464}, {425, 427, 463, 465}, {464, 428, 426, 466},
        {423, 425, 465, 467}, {466, 426, 424, 468}, {421, 423, 467, 469},
        {468, 424, 422, 470}, {419, 421, 469, 471}, {470, 422, 420, 472},
        {417, 419, 471, 473}, {472, 420, 418, 474}, {457, 455, 479, 477},
        {480, 456, 458, 478}, {477, 479, 481, 483}, {482, 480, 478, 484},
        {483, 481, 487, 485}, {488, 482, 484, 486}, {485, 487, 489, 491},
        {490, 488, 486, 492}, {463, 475, 485, 491}, {486, 476, 464, 492},
        {451, 483, 485, 475}, {486, 484, 452, 476}, {451, 461, 477, 483},
        {478, 462, 452, 484}, {457, 477, 461, 459}, {462, 478, 458, 460},
        {453, 473, 479, 455}, {480, 474, 454, 456}, {471, 481, 479, 473},
        {480, 482, 472, 474}, {469, 487, 481, 471}, {482, 488, 470, 472},
        {467, 489, 487, 469}, {488, 490, 468, 470}, {465, 491, 489, 467},
        {490, 492, 466, 468}, {391, 389, 503, 501}, {504, 390, 392, 502},
        {393, 391, 501, 499}, {502, 392, 394, 500}, {395, 393, 499, 497},
        {500, 394, 396, 498}, {397, 395, 497, 495}, {498, 396, 398, 496},
        {399, 397, 495, 493}, {496, 398, 400, 494}, {387, 399, 493, 505},
        {494, 400, 388, 506}, {493, 501, 503, 505}, {504, 502, 494, 506},
        {493, 495, 499, 501}, {500, 496, 494, 502}, {313, 381, 387, 505},
        {388, 382, 314, 506}, {313, 505, 503, 321}, {504, 506, 314, 322},
        {319, 321, 503, 389}, {504, 322, 320, 390}};
    auto shp = new shape{name};
    shp->pos = suzanne_pos;
    shp->quads = suzanne_quads;
    shp->quads.reserve(shp->quads.size() + suzanne_triangles.size());
    for (auto& t : suzanne_triangles) {
        shp->quads.push_back({t.x, t.y, t.z, t.z});
    }
    return shp;
}

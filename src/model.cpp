#include "scene.h"

#include "ext/yocto_math.h"
#include "ext/yocto_utils.h"
#include <algorithm>


#ifndef FALSE
	#define FALSE 0
#endif
#ifndef TRUE
	#define TRUE 1
#endif


#if !defined(DEBUGBREAK)
	#if defined(__clang__)
		#define DEBUGBREAK __builtin_trap()
	#elif defined(__GNUC__)
		#define DEBUGBREAK __asm__ __volatile__ ("int3")
	#elif defined(_MSC_VER)
		#define DEBUGBREAK __debugbreak()
	#endif
#endif

#if ENABLE_ASSERT
	#include <stdio.h>
	#define ASSERT(c) do{ if(!(c)) {fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); DEBUGBREAK;} } while(0)
#else
	#define ASSERT(c) do{ (void)sizeof(c); } while(0)
#endif


// converte un mesh in vertici sono condivis
// per ogni faccia, aggiunge i rispettivi vertici a mesh e aggiusta
// gli indici di faccia
// alla fine, fa l'update delle normali
// modifica la shape passata e ne ritorna il puntatore
shape*
facet_normals(shape* shp)
{
	std::vector<vec3f> poss  = std::vector<vec3f>();
	std::vector<vec3f> norms = std::vector<vec3f>();
	std::vector<vec3i> tris  = std::vector<vec3i>();
	
	for(int i = 0; i < shp->triangles.size(); ++i)
	{
		vec3f v0 = shp->pos[ shp->triangles[i]._[0] ];
		vec3f v1 = shp->pos[ shp->triangles[i]._[1] ];
		vec3f v2 = shp->pos[ shp->triangles[i]._[2] ];
		
		poss.push_back(v0);
		poss.push_back(v1);
		poss.push_back(v2);
		
		vec3f norm = triangle_normal(v0, v1, v2);
		norms.push_back(norm);
		norms.push_back(norm);
		norms.push_back(norm);
		
		tris.push_back(vec3i{3*i + 0, 3*i + 1, 3*i + 2});
	}
	
	std::swap(shp->pos,       poss);
	std::swap(shp->norm,      norms);
	std::swap(shp->triangles, tris);
	
	return shp;
}


// sposta ogni vertice i lungo la normale per una quantità disp_txt[i] * scale
// modifica la shape passata e ne ritorna il puntatore
shape*
displace(shape* shp, texture* disp_txt, float scale)
{
	for(int i = 0; i < shp->pos.size(); ++i)
	{
		vec3f displacement = eval_texture(disp_txt, shp->texcoord[i], false);
		shp->pos[i] += shp->norm[i] * displacement * scale;
	}
	compute_smooth_normals(shp);
	return shp;
}


// implementa catmull-clark subdivsion eseguendo la suddivisione level volte
// il primo step della suddivisione, che introduce i nuovi vertici e facce,
// e' implementato con tesselate()
// modifica la shape passata e ne ritorna il puntatore
shape*
catmull_clark(shape* shp, int level)
{
	for(int i = 0; i < level; ++i)
	{
		// step 1
		tesselate(shp);
		
		// step 2
		int vertices_count = shp->pos.size();
		vec3f* avg_v = new vec3f[vertices_count]();
		int*   avg_n = new   int[vertices_count]();
		
		for(int i = 0; i < shp->quads.size(); ++i)
		{
			vec4i face = shp->quads[i];
			vec3f centroid = 0.25f * (shp->pos[face.x] + shp->pos[face.y] + shp->pos[face.z] + shp->pos[face.w]);
			
			avg_v[face.x] += centroid;
			avg_v[face.y] += centroid;
			avg_v[face.z] += centroid;
			avg_v[face.w] += centroid;
			
			avg_n[face.x] += 1;
			avg_n[face.y] += 1;
			avg_n[face.z] += 1;
			avg_n[face.w] += 1;
		}
		
		for(int i = 0; i < vertices_count; ++i)
		{
			avg_v[i] *= 1.0f / avg_n[i];
		}
		
		// step 3
		for(int i = 0; i < vertices_count; ++i)
		{
			shp->pos[i] += (avg_v[i] - shp->pos[i]) * (4.0f / avg_n[i]);
		}
		
		delete[] avg_v;
		delete[] avg_n;
	}
	
	return shp;
}


// crea una quadrato nel piano (x,y) con dimensione [-r,r]x[-r,r]
// il piano ha usteps+1 in x and vsteps+1 in y
// le texture coordinates sono sono in [0,1]x[0,1]
// quindi ogni vertice ha pos=[-u*2+1,-v*2+1,1], norm=[0,0,1], texcoord=[u,v]
// name e' l'identificativo
// il mesh è fatto di quads con griglia definita da usteps x vsteps
// ogni quad e' definito come [i,j], [i+1,j], [i+1,j+1], [i,j+1]
shape*
make_quad(const std::string& name, int usteps, int vsteps, float r)
{
	shape* mesh = new shape{name};
	for(int y = 0; y <= vsteps; ++y)
	{
		for(int x = 0; x <= usteps; ++x)
		{
			float u = (float)x / usteps;
			float v = (float)y / vsteps;
			vec3f pos      = r * vec3f{u*2.0f - 1.0f, v*2.0f - 1.0f, 0.0f};
			vec3f norm     = {0, 0, 1};
			vec2f texcoord = {u, v};
			
			mesh->pos.push_back(pos);
			mesh->norm.push_back(norm);
			mesh->texcoord.push_back(texcoord);
		}
	}
	
	int stride = usteps+1;
	for(int y = 0; y < vsteps; ++y)
	{
		for(int x = 0; x < usteps; ++x)
		{
			int i0 =     y*stride + x;
			int i1 =     y*stride + x+1;
			int i2 = (y+1)*stride + x+1;
			int i3 = (y+1)*stride + x;
			vec4i indices = {i0, i1, i2, i3};
			mesh->quads.push_back(indices);
		}
	}
	return mesh;
}


// crea una sfera di raggio r in coordinate sferiche con (2 pi u, pi v)
// vedi la descrizione di make_quad()
// name e' l'identificativo
shape*
make_sphere(const std::string& name, int usteps, int vsteps, float r)
{
	shape* mesh = new shape{name};
	for(int y = 0; y <= vsteps; ++y)
	{
		for(int x = 0; x <= usteps; ++x)
		{
			float u = (float)x / usteps;
			float v = (float)y / vsteps;
			
			float phi   = 2.0f * pif * u;
			float theta = pif * v;
			
			vec3f pos      = r * vec3f{sinf(theta)*cosf(phi), sinf(theta)*sinf(phi), cosf(theta)};
			vec3f norm     = normalize(pos);
			vec2f texcoord = {u, 1.0f-v};
			
			mesh->pos.push_back(pos);
			mesh->norm.push_back(norm);
			mesh->texcoord.push_back(texcoord);
		}
	}
	
	int stride = usteps+1;
	for(int y = 0; y < vsteps; ++y)
	{
		for(int x = 0; x < usteps; ++x)
		{
			int i0 =     y*stride + x;
			int i1 =     y*stride + x+1;
			int i2 = (y+1)*stride + x+1;
			int i3 = (y+1)*stride + x;
			vec4i indices = {i0, i1, i2, i3};
			mesh->quads.push_back(indices);
		}
	}
	return mesh;
}


// crea una geosfera ottenuta tessellando la shape originale con tesselate()
// level volte e poi spostando i vertici risultanti sulla sfera di raggio r
// (i vertici sulla sfera unitaria sono normalize(p))
// name e' l'identificativo
shape*
make_geosphere(const std::string& name, int level, float r)
{
	const float X = 0.525731112119133606f;
	const float Z = 0.850650808352039932f;
	auto pos = std::vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z},
		{X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
		{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
	auto triangles = std::vector<vec3i>{{0, 1, 4}, {0, 4, 9}, {9, 4, 5},
		{4, 8, 5}, {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3}, {5, 3, 2},
		{2, 3, 7}, {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0}, {0, 6, 1},
		{6, 10, 1}, {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
	auto norm = pos;
	
	
	shape* mesh = new shape{name};
	mesh->pos = pos;
	mesh->norm = norm;
	mesh->triangles = triangles;
	for(int i = 0; i < level; ++i)
	{
		tesselate(mesh);
	}
	
	ASSERT(mesh->pos.size() == mesh->norm.size());
	for(int i = 0; i < mesh->pos.size(); ++i)
	{
		mesh->pos[i] = r * normalize(mesh->pos[i]);
		mesh->norm[i] = normalize(mesh->norm[i]);
	}
	
	return mesh;
}


// aggiunge un'instanza alla scena, assicurando che gli oggetti puntati
// dall'instanza (shp, mat, mat->kd_txt) sono anch'essi stati aggiunti alla
// scene almeno una e una sola volta
void
add_instance(scene* scn, const std::string& name, const frame3f& f, shape* shp, material* mat)
{
	if (!shp || !mat) return;
	
	if(std::find(scn->shapes.begin(), scn->shapes.end(), shp) == scn->shapes.end())
		scn->shapes.push_back(shp);
	
	if(std::find(scn->materials.begin(), scn->materials.end(), mat) == scn->materials.end())
		scn->materials.push_back(mat);
	
	if(mat->ke_txt && std::find(scn->textures.begin(), scn->textures.end(), mat->ke_txt) == scn->textures.end())
		scn->textures.push_back(mat->ke_txt);
	if(mat->kd_txt && std::find(scn->textures.begin(), scn->textures.end(), mat->kd_txt) == scn->textures.end())
		scn->textures.push_back(mat->kd_txt);
	if(mat->ks_txt && std::find(scn->textures.begin(), scn->textures.end(), mat->ks_txt) == scn->textures.end())
		scn->textures.push_back(mat->ks_txt);
	if(mat->rs_txt && std::find(scn->textures.begin(), scn->textures.end(), mat->rs_txt) == scn->textures.end())
		scn->textures.push_back(mat->rs_txt);
	if(mat->kr_txt && std::find(scn->textures.begin(), scn->textures.end(), mat->kr_txt) == scn->textures.end())
		scn->textures.push_back(mat->kr_txt);
	
	instance* inst = new instance{name, f, mat, shp};
	scn->instances.push_back(inst);
}


// crea an material con i parametrii specifici - name e' l'identificativo
material*
make_material(const std::string& name, const vec3f& kd, const std::string& kd_txt, const vec3f& ks = {0.04f, 0.04f, 0.04f}, float rs = 0.01f)
{
	material* mat = new material{name};
	mat->kd = kd;
	mat->ks = ks;
	mat->rs = rs;
	mat->kd_txt = new texture{kd_txt};;
	return mat;
}


// aggiunge num sfere di raggio r in un cerchio di raggio R orientato come
// specificato da f
void
add_sphere_instances(scene* scn, const frame3f& f, float R, float r, int num, material* mat)
{
	shape* mesh = make_sphere("sphere", 16, 16, r);
	for(int i = 0; i < num; ++i)
	{
		float angle = (2.0f * pif / num) * i;
		vec3f local_pos = R * vec3f{cosf(angle), sinf(angle), 0.0f};
		vec3f world_pos = transform_point(f, local_pos);
		frame3f frame = {f.x, f.y, f.z, world_pos};
		add_instance(scn, "sphere", frame, mesh, mat);
	}
}


scene* init_scene()
{
	auto scn = new scene();
	// add floor
	auto mat = new material{"floor"};
	mat->kd = {0.2f, 0.2f, 0.2f};
	mat->kd_txt = new texture{"grid.png"};
	scn->textures.push_back(mat->kd_txt);
	scn->materials.push_back(mat);
	auto shp = new shape{"floor"};
	shp->pos = {{-20, 0, -20}, {20, 0, -20}, {20, 0, 20}, {-20, 0, 20}};
	shp->norm = {{0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
	shp->texcoord = {{-10, -10}, {10, -10}, {10, 10}, {-10, 10}};
	shp->triangles = {{0, 1, 2}, {0, 2, 3}};
	scn->shapes.push_back(shp);
	scn->instances.push_back(new instance{"floor", identity_frame3f, mat, shp});
	// add light
	auto lshp = new shape{"light"};
	lshp->pos = {{1.4f, 8, 6}, {-1.4f, 8, 6}};
	lshp->points = {0, 1};
	scn->shapes.push_back(lshp);
	auto lmat = new material{"light"};
	lmat->ke = {100, 100, 100};
	scn->materials.push_back(lmat);
	scn->instances.push_back(new instance{"light", identity_frame3f, lmat, lshp});
	// add camera
	auto cam = new camera{"cam"};
	cam->frame = lookat_frame3f({0, 4, 10}, {0, 1, 0}, {0, 1, 0});
	cam->fovy = 15 * pif / 180.f;
	cam->aspect = 16.0f / 9.0f;
	cam->aperture = 0;
	cam->focus = length(vec3f{0, 4, 10} - vec3f{0, 1, 0});
	scn->cameras.push_back(cam);
	return scn;
}

int main(int argc, char** argv)
{
	// command line parsing
	auto parser =
		yu::cmdline::make_parser(argc, argv, "model", "creates simple scenes");
	auto sceneout = yu::cmdline::parse_opts(
		parser, "--output", "-o", "output scene", "out.obj");
	auto type = yu::cmdline::parse_args(
		parser, "type", "type fo scene to create", "empty", true);
	yu::cmdline::check_parser(parser);
	
	printf("creating scene %s\n", type.c_str());
	
	// create scene
	auto scn = init_scene();
	if (type == "empty") {
	} else if (type == "simple") {
		add_instance(scn, "quad", make_frame3_fromz({-1.25f, 1, 0}, {0, 0, 1}),
			make_quad("quad", 16, 16, 1),
			make_material("quad", {1, 1, 1}, "colored.png"));
		add_instance(scn, "sphere", make_frame3_fromz({1.25f, 1, 0}, {0, 0, 1}),
			make_sphere("sphere", 32, 16, 1),
			make_material("sphere", {1, 1, 1}, "colored.png"));
	} else if (type == "instances") {
		add_sphere_instances(scn,
			frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1.25f, 0}}, 1, 0.1, 16,
			make_material("obj", {1, 1, 1}, "colored.png"));
	} else if (type == "displace") {
		add_instance(scn, "quad1",
			make_frame3_fromz({-1.25f, 1, 0}, {0, 0, 1}),
			displace(make_quad("quad1", 64, 64, 1),
				make_grid_texture(256, 256),
				0.5),
			make_material("quad1", {1, 1, 1}, "colored.png"));
		add_instance(scn, "quad2",
			make_frame3_fromz({1.25f, 1, 0}, {0, 0, 1}),
			displace(make_quad("quad2", 64, 64, 1),
				make_bumpdimple_texture(256, 256),
				0.5),
			make_material("quad2", {1, 1, 1}, "colored.png"));
	} else if (type == "normals") {
		add_instance(scn, "smooth",
			make_frame3_fromz({-1.25f, 1, 0}, {0, 0, 1}),
			make_geosphere("smooth", 2, 1),
			make_material("smooth", {0.5f, 0.2f, 0.2f}, ""));
		add_instance(scn, "faceted",
			make_frame3_fromz({1.25f, 1, 0}, {0, 0, 1}),
			facet_normals(make_geosphere("faceted", 2, 1)),
			make_material("faceted", {0.2f, 0.5f, 0.2f}, ""));
	} else if (type == "subdiv") {
		add_instance(scn, "cube",
			make_frame3_fromzx({-1.25f, 1, 0}, {0, 0, 1}, {1, 0, 0}),
			catmull_clark(make_cube("cube"), 4),
			make_material("cube", {0.5f, 0.2f, 0.2f}, ""));
		add_instance(scn, "monkey",
			make_frame3_fromzx({1.25f, 1, 0}, {0, 0, 1}, {1, 0, 0}),
			catmull_clark(make_monkey("monkey"), 2),
			make_material("monkey", {0.2f, 0.5f, 0.2f}, ""));
	} else {
		throw std::runtime_error("bad scene type");
	}
	
	// save
	printf("saving scene %s\n", sceneout.c_str());
	save_scene(sceneout, scn);
	delete scn;
	return 0;
}

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <chrono>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <SFML/Graphics.hpp>
#include <tinyxml2/tinyxml2.h>

#include <selfup/ns_filesys.h>
#include <selfup/ns_helpers.h>

#define XASRT(b) do { if (! (b)) throw ErrExc(); } while (0)

// https://www.sfml-dev.org/tutorials/2.5/graphics-draw.php#drawing-from-threads
// https://www.sfml-dev.org/tutorials/2.5/window-window.php
//   window.setActive(false);  // window needs to deactivate so that an offthread may take over drawing
//   pollEvent/waitEvent       // events must be polled from thread, whom being the window's creator
// http://www.gamefromscratch.com/post/2015/10/26/SFML-CPP-Tutorial-Spritesheets-and-Animation.aspx
// https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/d#Path_commands
//   the m/M documentation is garbage, read it straight from spec (link below)
// https://svgwg.org/svg2-draft/paths.html#PathDataMovetoCommands
//   much superior m/M documentation
// https://stackoverflow.com/questions/1380371/what-are-the-most-widely-used-c-vector-matrix-math-linear-algebra-libraries-a#1452950
//   yay all vector libraries suck
// https://gamedevelopment.tutsplus.com/tutorials/collision-detection-using-the-separating-axis-theorem--gamedev-169
// (((!!isect)/2.0f)+0.5f)

using namespace tinyxml2;

template<typename T>
using sp = ::std::shared_ptr<T>;

typedef sf::Rect<float> Rectf;

class ErrExc : std::runtime_error
{
public:
	ErrExc() :
		std::runtime_error("")
	{}
};

class Img
{
public:
	Img(const char *fname, Rectf dim) :
		m_img(new sf::Image()),
		m_tex(new sf::Texture()),
		m_spr(),
		m_dim(dim)
	{
		if (! m_img->loadFromFile(fname))
			XASRT(0);
		if (! m_tex->loadFromImage(*m_img))
			XASRT(0);
		m_spr = sp<sf::Sprite>(new sf::Sprite(*m_tex));
		// SFML does something like
		//   0: start with identity transform
		//   1: transform scale around origin
		//   2: transform rotate around origin
		//   3: rawset translation
		//   4: rawset decrease translation by SCALED ORIGIN
		// this блять pants on head can be neutralized by
		// increasing translation to compensate.
		// we will increase by halfdim aka halfsize * scale.
		// since halfsize is origin, halfsize * scale is SCALED ORIGIN.
		const sf::Vector2u size(m_img->getSize());
		const sf::Vector2f halfsize(size.x / 2.0f, size.y / 2.0f);
		const sf::Vector2f halfdim(dim.width / 2.0f, dim.height / 2.0f);
		const sf::Vector2f scale(dim.width / size.x, dim.height / size.y);
		XASRT(fabsf(halfdim.x - halfsize.x * scale.x) < 0.001f && fabsf(halfdim.y - halfsize.y * scale.y) < 0.001f);
		m_spr->setOrigin(halfsize);
		m_spr->setScale(scale);
		m_spr->setPosition(sf::Vector2f(dim.left, dim.top) + halfdim);
	}

public:
	sp<sf::Image> m_img;
	sp<sf::Texture> m_tex;
	sp<sf::Sprite> m_spr;
	Rectf m_dim;
};

class Tri
{
public:
	Tri translatedBy(const sf::Vector2f &a)
	{
		Tri t;
		for (size_t i = 0; i < 3; i++)
			t.d[i] = d[i] + a;
		return t;
	}

	Rectf colRect()
	{
		const float mx = GS_MIN(d[0].x, GS_MIN(d[1].x, d[2].x));
		const float Mx = GS_MAX(d[0].x, GS_MAX(d[1].x, d[2].x));
		const float my = GS_MIN(d[0].y, GS_MIN(d[1].y, d[2].y));
		const float My = GS_MAX(d[0].y, GS_MAX(d[1].y, d[2].y));
		Rectf r(mx, my, Mx - mx, My - my);
		return r;
	}

public:
	sf::Vector2f d[3] = {};
};

class EntCol
{
public:
	virtual Tri* colTri(size_t a) = 0;
};

class E1 : public EntCol
{
public:
	E1(const sf::Vector2f &a, const sf::Vector2f &b, const sf::Vector2f &c) :
		m_t()
	{
		m_t.d[0] = a;
		m_t.d[1] = b;
		m_t.d[2] = c;
	}

	virtual Tri* colTri(size_t a) override
	{
		return a == 0 ? &m_t : nullptr;
	}

public:
	Tri m_t;
};

class QuadNode
{
public:
	// |0 1|
	// |2 3|
	sp<QuadNode> m_node[4] = {};
	std::map<sp<EntCol>, size_t> m_entry;
};

class QuadNode4
{
public:
	sp<QuadNode> m_d[4] = {};
};

class QuadTree
{
public:
	float m_bound;
	float m_min_bound;
	sp<QuadNode> m_root;
	std::map<sp<EntCol>, QuadNode4> m_ents;

	QuadTree(float bound, float min_bound) :
		m_bound(bound),
		m_min_bound(min_bound),
		m_root(new QuadNode())
	{}

	bool _boundContains(const Rectf &r)
	{
		if (r.left >= 0 && r.left + r.width < m_bound &&
			r.top >= 0 && r.top + r.height < m_bound)
			return true;
		return false;
	}

	QuadNode * _descendToEx(float inc_bound, float inc_x, float inc_y, const std::function<void(QuadNode *)> &cbextra, const std::function<sp<QuadNode>()> &cbnull)
	{
		float bas_bound = m_bound;
		float bas_top = 0, bas_left = 0;
		QuadNode *node = m_root.get();
		XASRT(bas_bound >= inc_bound);
		while (bas_bound != inc_bound) {
			bas_bound /= 2;
			const float midY = bas_top + bas_bound;
			const float midX = bas_left + bas_bound;
			const size_t nodeidx = (inc_y < midY ? 0 : 2) + (inc_x < midX ? 0 : 1);
			if (cbextra)
				cbextra(node);
			if (cbnull && ! node->m_node[nodeidx])
				node->m_node[nodeidx] = cbnull();
			// FIXME: probably should be an XASRT(0) ?
			if (! node->m_node[nodeidx])
				return nullptr;
			node = node->m_node[nodeidx].get();
			bas_top += inc_y < midY ? 0 : bas_bound;
			bas_left += inc_x < midX ? 0 : bas_bound;
		}
		return node;
	}

	QuadNode * _descendTo(float inc_bound, float inc_x, float inc_y)
	{
		return _descendToEx(inc_bound, inc_x, inc_y,
			nullptr,
			[&]() {
				return std::make_shared<QuadNode>();
			});
	}

	QuadNode * _descendToHarvestNocreate(float inc_bound, float inc_x, float inc_y, std::set<EntCol *> *cols)
	{
		return _descendToEx(inc_bound, inc_x, inc_y,
			[&](QuadNode *node) {
				for (auto it = node->m_entry.begin(); it != node->m_entry.end(); ++it)
					cols->insert(it->first.get());
			},
			nullptr);
	}

	void _floodHarvestNocreate(QuadNode *node, std::set<EntCol *> *cols)
	{
		if (! node)
			return;
		for (auto it = node->m_entry.begin(); it != node->m_entry.end(); ++it)
			cols->insert(it->first.get());
		for (size_t i = 0; i < 4; i++)
			_floodHarvestNocreate(node->m_node[i].get(), cols);
	}

	float _computeIncBound(const Rectf &r)
	{
		XASRT(_boundContains(r));
		const float maxside = GS_MAX(r.width, r.height);
		const float rank_ = ceil(log2f(maxside));
		const float inc_bound = exp2f(rank_);
		const float inc_bound_ = GS_MAX(inc_bound, m_min_bound);
		return inc_bound_;
	}

	void remove(const sp<EntCol> &ent)
	{
		auto it = m_ents.find(ent);
		if (it == m_ents.end())
			return;
		for (size_t i = 0; i < 4; i++) {
			if (! it->second.m_d[i])
				break;
			if (it->second.m_d[i]->m_entry.erase(ent) != 1)
				XASRT(0);
		}
	}

	void insert(const sp<EntCol> &ent, const Rectf &r)
	{
		const float inc_bound = _computeIncBound(r);
		// nodes may have duplicates
		QuadNode *nodes[4] = {
			_descendTo(inc_bound, r.left, r.top),            _descendTo(inc_bound, r.left + r.width, r.top),
			_descendTo(inc_bound, r.left, r.top + r.height), _descendTo(inc_bound, r.left + r.width, r.top + r.height),
		};
		for (size_t i = 0; i < 4; i++)
			nodes[i]->m_entry[ent] = 0;
	}

	void check(const Rectf &r, std::set<EntCol *> *o_cols)
	{
		std::set<EntCol *> cols;
		const float inc_bound = _computeIncBound(r);
		// nodes may have duplicates
		// FIXME: modify _descentToHarvestNocreate to take all four r corners at once
		//   and ensure no time is wasted harvesting same nodes twice
		QuadNode *nodes[4] = {
			_descendToHarvestNocreate(inc_bound, r.left, r.top, &cols), _descendToHarvestNocreate(inc_bound, r.left + r.width, r.top, &cols),
			_descendToHarvestNocreate(inc_bound, r.left, r.top + r.height, &cols), _descendToHarvestNocreate(inc_bound, r.left + r.width, r.top + r.height, &cols),
		};
		// FIXME: before calling _floodHarvestNocreate see to duplicates being filtered out
		//   ensuring no time is wasted harvesting same nodes twice
		for (size_t i = 0; i < 4; i++)
			_floodHarvestNocreate(nodes[i], &cols);
		*o_cols = std::move(cols);
	}
};

// FIXME:
#include <stuff.cpp>
// FIXME:

std::string killslash(const std::string &path)
{
	size_t off = path.size() - 1;
	while (off < path.size() && path[off] == '/' || path[off] == '\\')
		off--;
	return path.substr(0, off + 1);
}

std::string resource_path_find()
{
	const std::string &name = "misc.txt";
	const std::string &curdir = ns_filesys::current_executable_directory();
	const std::string &up0dir = killslash(curdir);
	const std::string &up1dir = killslash(ns_filesys::path_directory(up0dir));
	const std::string &up2dir = killslash(ns_filesys::path_directory(up1dir));
	std::fstream fs0(ns_filesys::path_append_abs_rel(up0dir, name), std::ios_base::in | std::ios_base::binary);
	std::fstream fs1(ns_filesys::path_append_abs_rel(up1dir, name), std::ios_base::in | std::ios_base::binary);
	std::fstream fs2(ns_filesys::path_append_abs_rel(up2dir, name), std::ios_base::in | std::ios_base::binary);
	if (fs0.good())
		return up0dir;
	else if (fs1.good())
		return up1dir;
	else if (fs2.good())
		return up2dir;
	else
		throw FilesysExc("resource path not found");
}

XMLElement * xml_elt_attr_val_find(const std::vector<XMLElement *> &v, const char *a, const char *x)
{
	for (size_t i = 0; i < v.size(); i++)
		if (v[i]->Attribute(a) && std::string(x).compare(v[i]->Attribute(a)) == 0)
			return v[i];
	return nullptr;
}

std::vector<XMLElement *> xml_child_sibling_find_all(XMLElement *e, const char *n)
{
	std::vector<XMLElement *> v;
	XMLElement *e2 = nullptr;
	if (!(e2 = e->FirstChildElement()))
		return v;
	if (std::string(n).compare(e2->Name()) == 0)
		v.push_back(e2);
	while ((e2 = e2->NextSiblingElement(n)))
		v.push_back(e2);
	return v;
}

std::vector<std::string> xml_paths_to_ds(const std::vector<XMLElement *> &paths)
{
	std::vector<std::string> v;
	for (size_t i = 0; i < paths.size(); i++) {
		const char *d = paths[i]->Attribute("d");
		XASRT(d);
		v.push_back(std::string(d));
	}
	return v;
}

sf::Vector2f xml_coordstr(const char *cs)
{
	float x, y;
	if (sscanf(cs, "%f,%f", &x, &y) != 2)
		XASRT(0);
	return sf::Vector2f(x, y);
}

Tri xml_d_to_tri(const char *d)
{
	std::vector<std::string> v;
	std::stringstream ss(d);
	std::string item;
	while (std::getline(ss, item, ' '))
		v.push_back(item);
	XASRT(ss.eof());
	XASRT(v.size() == 5 || v.size() == 6);
	XASRT(v.size() == 5); // FIXME:
	Tri t;
	if (v.size() == 5) {
		XASRT(v[4] == "z" || v[4] == "Z");
		XASRT(v[0] == "m" || v[0] == "M");
		t.d[0] = xml_coordstr(v[1].c_str());
		t.d[1] = xml_coordstr(v[2].c_str());
		t.d[2] = xml_coordstr(v[3].c_str());
		if (v[0] == "m") {
			t.d[1] += t.d[0];
			t.d[2] += t.d[1];
		}
	}
	return t;
}

std::vector<Tri> xml_ds_to_tris(const std::vector<std::string> &ds)
{
	std::vector<Tri> v;
	for (size_t i = 0; i < ds.size(); i++)
		v.push_back(xml_d_to_tri(ds[i].c_str()));
	return v;
}

sp<Img> xml_image_to_img(const XMLElement *im, const char *data_path)
{
	const char *href = im->Attribute("xlink:href");
	XASRT(href);
	std::string fname = ns_filesys::path_append_abs_rel(data_path, href);
	Rectf dim(im->FloatAttribute("x"), im->FloatAttribute("y"), im->FloatAttribute("width"), im->FloatAttribute("height"));
	return sp<Img>(new Img(fname.c_str(), dim));
}

std::vector<sp<Img> > xml_images_to_imgs(const std::vector<XMLElement *> ims, const char *data_path)
{
	std::vector<sp<Img> > v;
	for (size_t i = 0; i < ims.size(); i++)
		v.push_back(xml_image_to_img(ims[i], data_path));
	return v;
}

std::vector<sf::Vertex> verts_from_tris(const std::vector<Tri> &t)
{
	std::vector<sf::Vertex> v(t.size() * 3);
	for (size_t i = 0; i < t.size(); i++) {
		v[3 * i + 0].position = t[i].d[0];
		v[3 * i + 1].position = t[i].d[1];
		v[3 * i + 2].position = t[i].d[2];
		v[3 * i + 0].color = sf::Color(0, 255, 0);
		v[3 * i + 1].color = sf::Color(0, 255, 0);
		v[3 * i + 2].color = sf::Color(0, 255, 0);
	}
	return v;
}

float cross2(const Tri &points, const Tri &triangle)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/44269990#44269990
	auto pa = points.d[0];
	auto pb = points.d[1];
	auto pc = points.d[2];
	auto p0 = triangle.d[0];
	auto p1 = triangle.d[1];
	auto p2 = triangle.d[2];
	auto dXa = pa.x - p2.x;
	auto dYa = pa.y - p2.y;
	auto dXb = pb.x - p2.x;
	auto dYb = pb.y - p2.y;
	auto dXc = pc.x - p2.x;
	auto dYc = pc.y - p2.y;
	auto dX21 = p2.x - p1.x;
	auto dY12 = p1.y - p2.y;
	auto dX02 = p0.x - p2.x;
	auto dY20 = p2.y - p0.y;
	auto D = dY12 * dX02 - dX21 * dY20;
	auto sa = dY12 * dXa + dX21 * dYa;
	auto sb = dY12 * dXb + dX21 * dYb;
	auto sc = dY12 * dXc + dX21 * dYc;
	auto ta = dY20 * dXa + dX02 * dYa;
	auto tb = dY20 * dXb + dX02 * dYb;
	auto tc = dY20 * dXc + dX02 * dYc;
	if (D < 0) return ((sa >= 0 && sb >= 0 && sc >= 0) ||
		(ta >= 0 && tb >= 0 && tc >= 0) ||
		(sa+ta <= D && sb+tb <= D && sc+tc <= D));
	return ((sa <= 0 && sb <= 0 && sc <= 0) ||
		(ta <= 0 && tb <= 0 && tc <= 0) ||
		(sa+ta >= D && sb+tb >= D && sc+tc >= D));
}

bool triangles_intersect_4(const Tri &t0, const Tri &t1) {
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/44269990#44269990
	return !(cross2(t0,t1) ||
		cross2(t1,t0));
}

int main(int argc, char **argv)
{
	std::string resource_path = resource_path_find();
	std::string data_path = ns_filesys::path_append_abs_rel(resource_path, "data");

	XMLDocument doc;

	doc.LoadFile(ns_filesys::path_append_abs_rel(data_path, "test0.svg").c_str());
	XMLElement *doc_root = doc.RootElement();
	XASRT(doc_root);
	std::vector<XMLElement *> doc_g = xml_child_sibling_find_all(doc_root, "g");
	XMLElement *doc_g_collision_static = xml_elt_attr_val_find(doc_g, "id", "collision_static");
	XASRT(doc_g_collision_static);
	std::vector<XMLElement *> doc_g_collision_static_path = xml_child_sibling_find_all(doc_g_collision_static, "path");
	std::vector<std::string>  doc_g_collision_static_path_d = xml_paths_to_ds(doc_g_collision_static_path);
	std::vector<Tri>          doc_tris = xml_ds_to_tris(doc_g_collision_static_path_d);

	XMLElement *doc_g_background = xml_elt_attr_val_find(doc_g, "id", "background");
	XASRT(doc_g_background);
	std::vector<XMLElement *> doc_g_background_image = xml_child_sibling_find_all(doc_g_background, "image");

	std::vector<sp<Img> > imgs = xml_images_to_imgs(doc_g_background_image, data_path.c_str());

	std::vector<sf::Vertex> verts = verts_from_tris(doc_tris);

	sp<QuadTree> qt(new QuadTree(1024, 16));
	
	for (size_t i = 0; i < doc_tris.size(); i++) {
		sp<E1> e1(new E1(doc_tris[i].d[0], doc_tris[i].d[1], doc_tris[i].d[2]));
		qt->insert(e1, e1->m_t.colRect());
	}
	
	sf::RenderWindow window(sf::VideoMode(800, 600), "SF");

	window.setMouseCursorGrabbed(true);

	sf::VertexBuffer vb(sf::Triangles);
	if (! vb.create(verts.size()))
		XASRT(0);
	if (! vb.update(verts.data()))
		XASRT(0);

	Tri tb;
	tb.d[0] = sf::Vector2f(0, 0);
	tb.d[1] = sf::Vector2f(50, 0);
	tb.d[2] = sf::Vector2f(50, 50);
	sf::Vector2f tbpos;

	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				if (event.key.code == sf::Keyboard::Key::Escape)
					window.close();
				break;
			case sf::Event::MouseMoved:
				tbpos = sf::Vector2f((float)event.mouseMove.x, (float)event.mouseMove.y);
				break;
			}
		}
		window.clear(sf::Color(128, 128, 128));

		for (size_t i = 0; i < imgs.size(); i++)
			window.draw(*imgs[i]->m_spr);

		window.draw(vb);

		Tri q0 = tb.translatedBy(tbpos);
		
		std::vector<sf::Vertex> colverts;
		for (size_t i = 0; i < doc_tris.size(); i++) {
			std::set<EntCol *> cols1;
			qt->check(q0.colRect(), &cols1);
			for (auto it = cols1.begin(); it != cols1.end(); ++it) {
				const Tri &q1 = *(*it)->colTri(0);
				bool isect = triangles_intersect_4(q0, q1);
				if (isect)
					for (size_t j = 0; j < 3; j++)
						colverts.push_back(sf::Vertex(q1.d[j], sf::Color(255, 255, 0)));
			}
			//bool is2 = triangles_intersect_4(q0, doc_tris[i]);
			//if (is2)
			//	for (size_t j = 0; j < 3; j++)
			//		colverts.push_back(sf::Vertex(doc_tris[i].d[j], sf::Color(255, 255, 0)));
		}

		window.draw(colverts.data(), colverts.size(), sf::Triangles);

		std::vector<sf::Vertex> tvs;
		for (size_t i = 0; i < 3; i++)
			tvs.push_back(sf::Vertex(tb.d[i] + tbpos, sf::Color(0, 0, 255)));

		window.draw(tvs.data(), tvs.size(), sf::Triangles);

		window.display();
	}

	return EXIT_SUCCESS;
}

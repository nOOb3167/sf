#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <chrono>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
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

class Tri
{
public:
	sf::Vector2f d[3] = {};
};

class QuadNode
{
public:
	// |0 1|
	// |2 3|
	sp<QuadNode> m_node[4];
};

class QuadTree
{
public:
	float m_bound;
	float m_rank;
	sp<QuadNode> m_root;

	QuadTree(float bound) :
		m_bound(bound),
		m_rank(log2f(bound))
	{
		XASRT(m_rank == truncf(m_rank));
	}

	bool _bound_contains(const Rectf &r)
	{
		if (r.left >= 0 && r.left + r.width < m_bound &&
			r.top >= 0 && r.top + r.height < m_bound)
			return true;
		return false;
	}

	QuadNode * descendTo(size_t inc_bound, float inc_x, float inc_y)
	{
		size_t bas_bound = m_bound;
		float bas_top = 0, bas_left = 0;
		QuadNode *node = m_root.get();
		// descend thru nodes from bas_bound to inc_bound
		while (bas_bound != inc_bound) {
			bas_bound /= 2;
			const float midY = bas_top + bas_bound;
			const float midX = bas_left + bas_bound;
			node = node->m_node[(inc_y < midY ? 0 : 2) + (inc_x < midX ? 0 : 1)].get();
			bas_top += inc_y < midY ? 0 : bas_bound;
			bas_left += inc_x < midX ? 0 : bas_bound;
		}
		return node;
	}

	void insert(const Rectf &r)
	{
		XASRT(_bound_contains(r));
		const size_t maxside = GS_MAX(r.width, r.height);
		const size_t rank_ = truncf(log2f(maxside));
		XASRT(m_rank > rank_);
		const size_t rank = m_rank - rank_;

		//QuadNode *ul = m_root.get();
		QuadNode *parent = nullptr;
		for (size_t i = 0; i < rank; i++) {

		}

		//

		size_t inc_bound = exp2f(rank_);
		size_t bas_bound = m_bound;
		// descend thru nodes from bas_bound to inc_bound
		Rectf bas_r(0, 0, m_bound, m_bound);
		QuadNode *ul = m_root.get();
		QuadNode *ul_parent = nullptr;
		while (bas_bound != inc_bound) {
			const float midY = bas_r.top + bas_r.height / 2;
			const float midX = bas_r.left + bas_r.width / 2;
			ul_parent = ul;
			ul = ul->m_node[(r.top < midY ? 0 : 2) + (r.left < midX ? 0 : 1)].get();
			bas_r.top += r.top < midY ? 0 : bas_r.height / 2;
			bas_r.left += r.left < midX ? 0 : bas_r.width / 2;
			bas_r.height /= 2;
			bas_r.width /= 2;
			bas_bound /= 2;
		}

		//
		size_t inc_bound = exp2f(rank_);
		size_t bas_bound = m_bound;
		float bas_top = 0, bas_left = 0;
		QuadNode *ul = m_root.get();
		QuadNode *ul_parent = nullptr;
		// descend thru nodes from bas_bound to inc_bound
		while (bas_bound != inc_bound) {
			bas_bound /= 2;
			const float midY = bas_top + bas_bound;
			const float midX = bas_left + bas_bound;
			ul_parent = ul;
			ul = ul->m_node[(r.top < midY ? 0 : 2) + (r.left < midX ? 0 : 1)].get();
			bas_top += r.top < midY ? 0 : bas_bound;
			bas_left += r.left < midX ? 0 : bas_bound;
		}
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

	XMLDocument doc;

	doc.LoadFile(ns_filesys::path_append_abs_rel(resource_path, "data/test0.svg").c_str());
	XMLElement *doc_root = doc.RootElement();
	XASRT(doc_root);
	std::vector<XMLElement *> doc_g = xml_child_sibling_find_all(doc_root, "g");
	std::vector<XMLElement *> doc_g_path = xml_child_sibling_find_all(doc_g.at(0), "path");
	std::vector<std::string>  doc_d = xml_paths_to_ds(doc_g_path);
	std::vector<Tri>          doc_tris = xml_ds_to_tris(doc_d);

	std::vector<sf::Vertex> verts = verts_from_tris(doc_tris);

	sf::RenderWindow window(sf::VideoMode(800, 600), "SF");

	window.setMouseCursorGrabbed(true);

	sf::RectangleShape shape(sf::Vector2f(16, 16));
	shape.setFillColor(sf::Color(255, 0, 0));

	sf::VertexBuffer vb(sf::Triangles);
	if (! vb.create(verts.size()))
		XASRT(0);
	if (! vb.update(verts.data()))
		XASRT(0);

	sf::Vertex tb[3];
	tb[0].position = sf::Vector2f(0, 0);
	tb[1].position = sf::Vector2f(50, 0);
	tb[2].position = sf::Vector2f(50, 50);
	sf::Vertex tvs[6];
	tvs[0].position = sf::Vector2f(400, 400);
	tvs[1].position = sf::Vector2f(500, 350);
	tvs[2].position = sf::Vector2f(450, 450);
	tvs[3].position = tb[0].position;
	tvs[4].position = tb[1].position;
	tvs[5].position = tb[2].position;
	for (size_t i = 0; i < 6; i++)
		tvs[i].color = sf::Color(0, 0, 255);

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
				sf::Vector2f mcenter((float)event.mouseMove.x, (float)event.mouseMove.y);
				tvs[3].position = tb[0].position + mcenter;
				tvs[4].position = tb[1].position + mcenter;
				tvs[5].position = tb[2].position + mcenter;
				break;
			}
		}
		window.clear(sf::Color(128, 128, 128));

		window.draw(shape, sf::RenderStates(sf::Transform().translate(16, 16)));
		window.draw(vb);

		Tri q0;
		q0.d[0] = tvs[0].position;
		q0.d[1] = tvs[1].position;
		q0.d[2] = tvs[2].position;
		Tri q1;
		q1.d[0] = tvs[3].position;
		q1.d[1] = tvs[4].position;
		q1.d[2] = tvs[5].position;
		bool isect = triangles_intersect_4(q0, q1);
		for (size_t i = 0; i < 6; i++)
			tvs[i].color = sf::Color(0, 0, 255 * (((!!isect)/2.0f)+0.5f));

		window.draw(tvs, 6, sf::Triangles);

		window.display();
	}

	return EXIT_SUCCESS;
}
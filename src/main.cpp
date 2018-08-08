#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <chrono>
#include <fstream>
#include <functional>
#include <limits>
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

#define DF_ENTCOL_MAX_TRIS 20
#define DF_ENTCOL_MAX_MOVE 20.0f
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

/* probably not required anyway just remove if it fails on any platform */
static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");

using namespace tinyxml2;

template<typename T>
using sp = ::std::shared_ptr<T>;

const auto rel = ::ns_filesys::path_append_abs_rel;

typedef sf::Rect<float> Rectf;

class Tri;

bool triangles_intersect_4(const Tri &t0, const Tri &t1);
bool triangles_intersect_4_offset(const sf::Vector2f &t0_offset, const Tri &t0, const Tri &t1);

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

	static std::pair<Tri, Tri> fromRect(const Rectf &dim)
	{
		std::pair<Tri, Tri> t2;
		Tri &t0 = t2.first, &t1 = t2.second;
		t0.d[0] = { dim.left, dim.top };
		t0.d[1] = { dim.left + dim.width, dim.top };
		t0.d[2] = { dim.left + dim.width, dim.top + dim.height };
		t1.d[0] = { dim.left, dim.top };
		t1.d[1] = { dim.left, dim.top + dim.height };
		t1.d[2] = { dim.left + dim.width, dim.top + dim.height };
		return t2;
	}

public:
	sf::Vector2f d[3] = {};
};

class Img
{
public:
	Img(const char *fname, Rectf dim) :
		m_img(new sf::Image()),
		m_tex(new sf::Texture()),
		m_spr(),
		m_dim(dim),
		m_dim_tris(2)
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
		const sf::Vector2f halfdim(m_dim.width / 2.0f, m_dim.height / 2.0f);
		const sf::Vector2f scale(m_dim.width / size.x, m_dim.height / size.y);
		XASRT(fabsf(halfdim.x - halfsize.x * scale.x) < 0.001f && fabsf(halfdim.y - halfsize.y * scale.y) < 0.001f);
		m_spr->setOrigin(halfsize);
		m_spr->setScale(scale);
		m_spr->setPosition(sf::Vector2f(m_dim.left, m_dim.top) + halfdim);
		_refreshDimTris();
	}

	void _refreshDimTris()
	{
		auto t2 = Tri::fromRect(m_dim);
		m_dim_tris[0] = t2.first;
		m_dim_tris[1] = t2.second;
	}

	void setPosition(float x, float y)
	{
		const sf::Vector2f halfdim(m_dim.width / 2.0f, m_dim.height / 2.0f);
		m_spr->setPosition(sf::Vector2f(x, y) + halfdim);
		m_dim.left = x;
		m_dim.top = y;
		_refreshDimTris();
	}

	Rectf & getDim()
	{
		return m_dim;
	}

public:
	sp<sf::Image> m_img;
	sp<sf::Texture> m_tex;
	sp<sf::Sprite> m_spr;
	Rectf m_dim;
	std::vector<Tri> m_dim_tris;
};

class EntCol
{
public:
	virtual const std::vector<Tri> & colTri() const = 0;

	Rectf colRect() const
	{
		float mx = std::numeric_limits<float>::infinity();
		float Mx = -std::numeric_limits<float>::infinity();
		float my = std::numeric_limits<float>::infinity();
		float My = -std::numeric_limits<float>::infinity();
		for (auto it = colTri().begin(); it != colTri().end(); ++it) {
			mx = GS_MIN(mx, GS_MIN(it->d[0].x, GS_MIN(it->d[1].x, it->d[2].x)));
			Mx = GS_MAX(Mx, GS_MAX(it->d[0].x, GS_MAX(it->d[1].x, it->d[2].x)));
			my = GS_MIN(my, GS_MIN(it->d[0].y, GS_MIN(it->d[1].y, it->d[2].y)));
			My = GS_MAX(My, GS_MAX(it->d[0].y, GS_MAX(it->d[1].y, it->d[2].y)));
		}
		Rectf r(mx, my, Mx - mx, My - my);
		return r;
	}
};

typedef ::std::map<sp<EntCol>, size_t> ent_map_t;

class ImgCurve
{
public:
	ImgCurve(const char *fname) :
		m_img(new sf::Image()),
		m_size(),
		m_pts()
	{
		if (! m_img->loadFromFile(fname))
			XASRT(0);
		m_size = sf::Vector2f(m_img->getSize().x, m_img->getSize().y);
		m_pts = std::vector<float>(m_size.x);
		float ym1 = m_size.y - 1;
		for (size_t x = 0; x < m_size.x; x++) {
			float num_transparent = 0;
			for (size_t y = 0; y < m_size.y; y++)
				if (m_img->getPixel(x, ym1 - y).a == 0)
					num_transparent++;
				else
					break;
			m_pts[x] = num_transparent;
		}
		m_pts_inc = std::vector<float>(m_size.x);
		float cur_lvl = 0;
		for (size_t x = 0; x < m_size.x; x++) {
			m_pts_inc[x] = m_pts[x] - cur_lvl;
			cur_lvl = m_pts[x];
		}
	}

public:
	sp<sf::Image> m_img;
	sf::Vector2f m_size;
	std::vector<float> m_pts;
	std::vector<float> m_pts_inc;
};

class E1 : public EntCol
{
public:
	E1(const sf::Vector2f &a, const sf::Vector2f &b, const sf::Vector2f &c) :
		m_t(1)
	{
		m_t[0].d[0] = a;
		m_t[0].d[1] = b;
		m_t[0].d[2] = c;
	}

	virtual const std::vector<Tri> & colTri() const override
	{
		return m_t;
	}

public:
	std::vector<Tri> m_t;
};

class E2 : public EntCol
{
public:
	E2(const std::string &data_path, const Rectf &dim) :
		m_img(new Img(rel(data_path, "e2_0.png").c_str(), dim))
	{}

	virtual const std::vector<Tri> & colTri() const override
	{
		return m_img->m_dim_tris;
	}

public:
	sp<Img> m_img;
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
	QuadNode * m_d[4] = {};
};

class QuadTree
{
public:
	float m_bound;
	float m_min_bound;
	sp<QuadNode> m_root;
	std::map<sp<EntCol>, QuadNode4> m_ents;

	sf::Transform m_cw;
	sf::Transform m_ccw;

	QuadTree(float bound, float min_bound) :
		m_bound(bound),
		m_min_bound(min_bound),
		m_root(new QuadNode()),
		m_cw(sf::Transform().rotate(15)),
		m_ccw(sf::Transform().rotate(-15))
	{}

	template <typename InputIter>
	static QuadTree * createFromEntCol(float bound, float min_bound, InputIter &esb, InputIter &ese)
	{
		QuadTree * qt = new QuadTree(bound, min_bound);
		for (/*dummy*/; esb != ese; ++esb)
			qt->insert(esb->first, esb->first->colRect());
		return qt;
	}

	bool _boundContains(const Rectf &r) const
	{
		if (r.left >= 0 && r.left + r.width < m_bound &&
			r.top >= 0 && r.top + r.height < m_bound)
			return true;
		return false;
	}

	QuadNode * _descendTo(float inc_bound, float inc_x, float inc_y) const
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
			//
			if (! node->m_node[nodeidx])
				node->m_node[nodeidx] = std::make_shared<QuadNode>();
			//
			node = node->m_node[nodeidx].get();
			bas_top += inc_y < midY ? 0 : bas_bound;
			bas_left += inc_x < midX ? 0 : bas_bound;
		}
		return node;
	}

	QuadNode * _descendToHarvestNocreate(float inc_bound, float inc_x, float inc_y, std::set<EntCol *> *cols) const
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
			//
			for (auto it = node->m_entry.begin(); it != node->m_entry.end(); ++it)
				cols->insert(it->first.get());
			if (! node->m_node[nodeidx])
				return nullptr;
			//
			node = node->m_node[nodeidx].get();
			bas_top += inc_y < midY ? 0 : bas_bound;
			bas_left += inc_x < midX ? 0 : bas_bound;
		}
		return node;
	}

	void _floodHarvestNocreate(QuadNode *node, std::set<EntCol *> *cols) const
	{
		if (! node)
			return;
		for (auto it = node->m_entry.begin(); it != node->m_entry.end(); ++it)
			cols->insert(it->first.get());
		for (size_t i = 0; i < 4; i++)
			_floodHarvestNocreate(node->m_node[i].get(), cols);
	}

	float _computeIncBound(const Rectf &r) const
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
		XASRT(m_ents.find(ent) == m_ents.end());
		m_ents[ent] = QuadNode4();
		m_ents[ent].m_d[0] = nodes[0];
		m_ents[ent].m_d[1] = nodes[1];
		m_ents[ent].m_d[2] = nodes[2];
		m_ents[ent].m_d[3] = nodes[3];
	}

	void check(const Rectf &r, std::set<EntCol *> *o_cols) const
	{
		const float inc_bound = _computeIncBound(r);
		// nodes may have duplicates
		// FIXME: modify _descentToHarvestNocreate to take all four r corners at once
		//   and ensure no time is wasted harvesting same nodes twice
		QuadNode *nodes[4] = {
			_descendToHarvestNocreate(inc_bound, r.left, r.top, o_cols), _descendToHarvestNocreate(inc_bound, r.left + r.width, r.top, o_cols),
			_descendToHarvestNocreate(inc_bound, r.left, r.top + r.height, o_cols), _descendToHarvestNocreate(inc_bound, r.left + r.width, r.top + r.height, o_cols),
		};
		// FIXME: before calling _floodHarvestNocreate see to duplicates being filtered out
		//   ensuring no time is wasted harvesting same nodes twice
		for (size_t i = 0; i < 4; i++)
			_floodHarvestNocreate(nodes[i], o_cols);
	}

	bool checkFreeSpot(const sf::Vector2f &offset_ec, const std::vector<Tri> &t_ec, const Rectf &r_ec, const std::set<EntCol *> &cols) const
	{
		const Rectf r(r_ec.left + offset_ec.x, r_ec.top + offset_ec.y, r_ec.width, r_ec.height);
		for (auto it = cols.begin(); it != cols.end(); ++it) {
			const Rectf &rr = (*it)->colRect();
			// rect-level collision
			if (! r.intersects(rr))
				continue;
			// tris-level collision
			for (auto it2 = t_ec.begin(); it2 != t_ec.end(); ++it2)
				for (auto it3 = (*it)->colTri().begin(); it3 != (*it)->colTri().end(); ++it3)
					if (triangles_intersect_4_offset(offset_ec, *it2, *it3))
						return false;
		}
		return true;
	}

	sf::Vector2f checkSimulate(const EntCol &ec, const sf::Vector2f &move, const sf::Vector2f &extramove) const
	{
		// FIXME: imagine ec just above an obstacle and moving down by 5.
		//   ends up embedded by 4 inside obstacle.
		//   cw,ccw must be practically 180deg reverse to backtrack out.
		//   recommend adding a scale factor to cw,ccw.
		const std::vector<Tri> &t_ec = ec.colTri();
		const Rectf &r_ec = ec.colRect();
		// movements max_move limited per-axis not by vector magnitude
		const sf::Vector2f cur(r_ec.left, r_ec.top);
		const sf::Vector2f &mv = move + extramove;
		const sf::Vector2f &twicemv = mv * 2.0f;
		const sf::Vector2f &nxt = cur + mv;
		XASRT(fabsf(mv.x) < DF_ENTCOL_MAX_MOVE && fabsf(mv.y) < DF_ENTCOL_MAX_MOVE);
		const Rectf r_nxt_expanded(nxt.x - mv.x, nxt.y - mv.y, r_ec.width + twicemv.x, r_ec.height + twicemv.y);
		std::set<EntCol *> cols;
		// FIXME: (perf) possibly check vs r_nxt and go for r_nxt_expanded only on a hit
		check(r_nxt_expanded, &cols);
		if (cols.empty())
			return nxt;
		if (checkFreeSpot(mv, t_ec, r_ec, cols))
			return nxt;
		sf::Vector2f cw = mv;
		sf::Vector2f ccw = mv;
		for (size_t i = 0; i < 3; i++) {
			cw = m_cw.transformPoint(cw);
			ccw = m_ccw.transformPoint(ccw);
			if (checkFreeSpot(cw, t_ec, r_ec, cols))
				return cur + cw;
			if (checkFreeSpot(ccw, t_ec, r_ec, cols))
				return cur + ccw;
		}
		// no free spot - stay still
		return cur;
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
	std::fstream fs0(rel(up0dir, name), std::ios_base::in | std::ios_base::binary);
	std::fstream fs1(rel(up1dir, name), std::ios_base::in | std::ios_base::binary);
	std::fstream fs2(rel(up2dir, name), std::ios_base::in | std::ios_base::binary);
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
	Rectf dim(im->FloatAttribute("x"), im->FloatAttribute("y"), im->FloatAttribute("width"), im->FloatAttribute("height"));
	return sp<Img>(new Img(rel(data_path, href).c_str(), dim));
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

int main(int argc, char **argv)
{
	std::string resource_path = resource_path_find();
	std::string data_path = rel(resource_path, "data");

	XMLDocument doc;

	doc.LoadFile(rel(data_path, "test0.svg").c_str());
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

	ent_map_t es;

	for (size_t i = 0; i < doc_tris.size(); i++)
		es[sp<EntCol>(new E1(doc_tris[i].d[0], doc_tris[i].d[1], doc_tris[i].d[2]))] = 0;

	sp<const QuadTree> qt_static(QuadTree::createFromEntCol(1024, 16, es.begin(), es.end()));
	
	sp<E2> e2(new E2(data_path, Rectf(0, 0, 50, 50)));
	sp<E2> e2fall(new E2(data_path, Rectf(400, 0, 50, 50)));

	sp<ImgCurve> cur0(new ImgCurve(rel(data_path, "cur0.png").c_str()));

	sf::RenderWindow window(sf::VideoMode(800, 600), "SF");

	window.setMouseCursorGrabbed(true);
	window.setFramerateLimit(60);

	sf::VertexBuffer vb(sf::Triangles);
	if (! vb.create(verts.size()))
		XASRT(0);
	if (! vb.update(verts.data()))
		XASRT(0);

	Tri tb;
	tb.d[0] = sf::Vector2f(0, 0);
	tb.d[1] = sf::Vector2f(50, 0);
	tb.d[2] = sf::Vector2f(50, 50);

	sf::Vector2f curgravity(0, 5);
	size_t jumping = -1;

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
				e2->m_img->setPosition((float)event.mouseMove.x, (float)event.mouseMove.y);
				break;
			}
		}
		window.clear(sf::Color(128, 128, 128));

		for (size_t i = 0; i < imgs.size(); i++)
			window.draw(*imgs[i]->m_spr);

		window.draw(vb);

		sf::Vector2f wantmove;
		bool wantjump = false;

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A))
			wantmove.x -= 5;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D))
			wantmove.x += 5;
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space))
			wantjump = true;

		if (wantjump && jumping == -1)
			jumping = 0;
		if (jumping == cur0->m_size.x) {
			jumping = -1;
			curgravity = sf::Vector2f(0, 5);
		}
		if (jumping >= 0 && jumping < cur0->m_size.x) {
			wantmove.y -= cur0->m_pts_inc[jumping];
			jumping++;
			curgravity = sf::Vector2f(0, 0);
		}

		{
			const sf::Vector2f &newpos = qt_static->checkSimulate(*e2fall, wantmove, curgravity);
			e2fall->m_img->setPosition(newpos.x, newpos.y);
		}

		std::vector<sf::Vertex> colverts;
		{
			std::set<EntCol *> cols1;
			qt_static->check(e2->colRect(), &cols1);
			for (auto it = cols1.begin(); it != cols1.end(); ++it) {
				// FIXME: obsolete colTri 0 code
				const Tri &q1 = (*it)->colTri()[0];
				bool isect = triangles_intersect_4(e2->colTri()[0], q1);
				if (isect)
					for (size_t j = 0; j < 3; j++)
						colverts.push_back(sf::Vertex(q1.d[j], sf::Color(255, 255, 0)));
			}
		}

		window.draw(colverts.data(), colverts.size(), sf::Triangles);

		window.draw(*e2->m_img->m_spr);
		window.draw(*e2fall->m_img->m_spr);

		window.display();
	}

	return EXIT_SUCCESS;
}

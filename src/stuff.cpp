bool sameside_broken(const sf::Vector2f &p1_, const sf::Vector2f &p2_, const sf::Vector2f &a_, const sf::Vector2f &b_)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/38550586#38550586
	const Eigen::Vector3f p1(p1_.x, p1_.y, 0);
	const Eigen::Vector3f p2(p2_.x, p2_.y, 0);
	const Eigen::Vector3f a(a_.x, a_.y, 0);
	const Eigen::Vector3f b(b_.x, b_.y, 0);
	const Eigen::Vector3f bma = b - a;
	const Eigen::Vector3f cp1 = bma.cross(p1 - a);
	const Eigen::Vector3f cp2 = bma.cross(p2 - a);
	return cp1.dot(cp2) >= 0;
}

bool triangle_intersection_broken(const Tri &P1, const Tri &P2)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/38550586#38550586
	for (size_t i = 0; i < 3; i++)
		if (! sameside_broken(P1.d[i], P2.d[0], P1.d[(i+1)%3], P1.d[(i+2)%3])
			&& sameside_broken(P2.d[0], P2.d[1], P1.d[(i+1)%3], P1.d[(i+2)%3])
			&& sameside_broken(P2.d[1], P2.d[2], P1.d[(i+1)%3], P1.d[(i+2)%3]))
		{
			return false;
		}
	return true;
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

bool triangles_intersect_4(const Tri &t0, const Tri &t1)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/44269990#44269990
	return !(cross2(t0,t1) ||
		cross2(t1,t0));
}

float cross2_offset(const sf::Vector2f &t0_offset, const Tri &points, const Tri &triangle)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/44269990#44269990
	auto pa = points.d[0] + t0_offset;
	auto pb = points.d[1] + t0_offset;
	auto pc = points.d[2] + t0_offset;
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

bool triangles_intersect_4_offset(const sf::Vector2f &t0_offset, const Tri &t0, const Tri &t1)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/44269990#44269990
	return !(cross2_offset(t0_offset, t0,t1) ||
		cross2_offset(t0_offset, t1,t0));
}

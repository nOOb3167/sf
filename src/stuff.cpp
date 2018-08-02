bool sameside(const sf::Vector2f &p1_, const sf::Vector2f &p2_, const sf::Vector2f &a_, const sf::Vector2f &b_)
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

bool triangle_intersection(const Tri &P1, const Tri &P2)
{
	// https://stackoverflow.com/questions/2778240/detection-of-triangle-collision-in-2d-space/38550586#38550586
	for (size_t i = 0; i < 3; i++)
		if (! sameside(P1.d[i], P2.d[0], P1.d[(i+1)%3], P1.d[(i+2)%3])
			&& sameside(P2.d[0], P2.d[1], P1.d[(i+1)%3], P1.d[(i+2)%3])
			&& sameside(P2.d[1], P2.d[2], P1.d[(i+1)%3], P1.d[(i+2)%3]))
		{
			return false;
		}
	return true;
}

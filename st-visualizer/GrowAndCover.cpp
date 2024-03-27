
#include "GrowAndCover.h"

#include <map>

using std::map;
using std::pair;
using std::vector;

// Construct a 2D counterclockwise rotation from the angle in radian.
Eigen::Rotation2Df rotM(float a)
{
	return Eigen::Rotation2D(a);
}

// ReSharper disable once CppInconsistentNaming
Eigen::Rotation2D hexM = rotM(pi / 3);

// Pull out all the points which lie on the grid. Only works in hex space.
// v2 is pi/3 radians from v1 counterclockwise with the same magnitude.
std::pair<std::vector<int>, Eigen::Matrix2Xi> getInliers(const Eigen::Matrix2Xf &pts, const Eigen::Vector2f &origin,
														 const Eigen::Vector2f &v1)
{
	const Eigen::Vector2f v2 = hexM * v1;
	Eigen::Matrix2Xi int_coords = roundPtsToCoords(pts, origin, v1, v2);
	const float errorMargin = HEX_ROUNDING_ERROR * v1.norm();

	// Check which indices are within the appropriate bounds
	std::vector<int> indices;
	indices.reserve(pts.cols());
	for (int i = 0; i < pts.cols(); i++)
	{
		const auto &pt = pts.col(i);
		Eigen::Vector2f returnPt = getPoint(int_coords.col(i).cast<float>(), origin, v1, v2);

		const float realMargin = (pt - returnPt).norm();
		if (realMargin < errorMargin)
		{
			indices.push_back(i);
		}
	}

	Eigen::Matrix2Xi revisedCoords = Eigen::Matrix2Xi::Zero(2, static_cast<int>(indices.size()));
	for (size_t i = 0; i < indices.size(); i++)
	{
		revisedCoords.col(i) = int_coords.col(indices[i]);
	}

	return {indices, revisedCoords};
}

Eigen::Matrix2Xi roundPtsToCoords(const Eigen::Matrix2Xf &pts, const Eigen::Vector2f &origin, const Eigen::Vector2f &v1, const Eigen::Vector2f &v2)
{
	Eigen::Matrix2Xi coords = Eigen::Matrix2Xi::Zero(2, pts.cols());
	for (long long i = 0; i < pts.cols(); i++)
	{
		coords.col(i) = roundPtToCoord(pts.col(i), origin, v1, v2);
	}
	return coords;
}

// returns {best origin, best v1},resulting inliers
pair<pair<Eigen::Vector2f, Eigen::Vector2f>, pair<vector<int>, Eigen::Matrix2Xi>> initGridInliers(
	const Eigen::Matrix2Xf &pts, const int &num)
{
	if (pts.cols() < 2)
		throw "Need more points to test";

	std::pair<std::vector<int>, Eigen::Matrix2Xi> best_inliers;
	Eigen::Vector2f best_origin({0, 0});
	Eigen::Vector2f best_v1({1, 0});

	for (int i = 0; i < num; i++)
	{
		// Pick a random point
		const long long origin_index = rand() % pts.cols();
		Eigen::Vector2f rand_ori = pts.col(origin_index);

		// Get the next closest point
		Eigen::Matrix2Xf adjusted_points = pts.colwise() - rand_ori;
		int selected_point = origin_index > 0 ? 0 : 1;

		float min_norm = adjusted_points.col(selected_point).norm();
		for (int j = 0; j < adjusted_points.cols(); j++)
		{
			const float test_norm = adjusted_points.col(j).norm();
			if (j != origin_index && test_norm < min_norm)
			{
				selected_point = j;
				min_norm = test_norm;
			}
		}

		const Eigen::Vector2f v1 = adjusted_points.col(selected_point);

		// See how many inliers there are with the resulting grid
		auto inliers = getInliers(pts, rand_ori, v1);

		// If it's a better fit, save it
		if (inliers.first.size() > best_inliers.first.size())
		{
			best_inliers = inliers;
			best_origin = rand_ori;
			best_v1 = v1;
		}
	}

	return std::pair(std::pair(best_origin, best_v1), best_inliers);
}

std::pair<Eigen::Vector2f, Eigen::Vector2f> getGrid(const Eigen::Matrix2Xf &pts, const std::vector<int> &indices, const Eigen::Matrix2Xi &intCoords)
{
	Eigen::Matrix4f ata = Eigen::Matrix4f::Zero();
	Eigen::Vector4f atb = Eigen::Vector4f::Zero();
	const Eigen::Matrix2f identity2 = Eigen::Matrix2f::Identity();

	for (int i = 0; i < intCoords.cols(); i++)
	{
		// For each coordinate
		const Eigen::Block<const Eigen::Matrix<int, 2, -1, 0, 2, -1>, 2, 1, 1> coord = intCoords.col(i);

		Eigen::Matrix<float, 2, 4> a;
		a << coord(0) * identity2 + coord(1) * hexM.toRotationMatrix(), identity2;
		ata += a.transpose() * a;
		atb += a.transpose() * pts.col(indices[i]);
	}

	Eigen::Vector4f result = ata.inverse() * atb;
	Eigen::Vector2f origin({result(2), result(3)});
	Eigen::Vector2f v1({result(0), result(1)});

	return {origin, v1};
}

// Refine Grid
std::pair<Eigen::Vector2f, Eigen::Vector2f> refineGrid(const Eigen::Matrix2Xf &pts,
													   const std::pair<Eigen::Vector2f, Eigen::Vector2f> &grid,
													   const std::pair<vector<int>, Eigen::Matrix2Xi> &inliers)
{
	pair<Eigen::Vector2f, Eigen::Vector2f> new_grid = grid;
	pair<vector<int>, Eigen::Matrix2Xi> new_inliers = inliers;
	size_t num = inliers.first.size();

	const std::function std_dev([](const Eigen::Matrix2Xf &_pts, const std::pair<Eigen::Vector2f, Eigen::Vector2f> &_grid)
								{
		Eigen::Matrix2Xf delta = _pts;
		const auto& origin = _grid.first;
		const auto& v1 = _grid.second;
		const auto v2 = hexM * _grid.second;
		for(int i = 0; i < _pts.cols(); i++)
		{
			delta.col(i) -= getPoint(roundPtToCoord(_pts.col(i),origin, v1, v2).cast<float>(), origin, v1, v2);
		}
		return delta.colwise().squaredNorm().rowwise().sum().col(0)(0); });

	float s = std_dev(pts, new_grid);
	while (static_cast<long long>(num) < pts.cols())
	{

		// Get the best origin and v1 from the grid
		new_grid = getGrid(pts, new_inliers.first, new_inliers.second);

		// Then get the inliers and grid that results from the origin and v1
		new_inliers = getInliers(pts, new_grid.first, new_grid.second);
		// If there is an improvement, go again
		const float new_s = std_dev(pts, new_grid);

		if (new_inliers.first.size() > num)
		{
			num = new_inliers.first.size();
			s = new_s;
		}
		else
		{
			break;
		}
	}

	new_grid = getGrid(pts, new_inliers.first, new_inliers.second);
	return new_grid;
}

pair<pair<Eigen::Vector2f, Eigen::Vector2f>, Eigen::Matrix2Xi> getGridAndCoords(const Eigen::Matrix2Xf &pts, const int &num)
{
	const pair<pair<Eigen::Vector2f, Eigen::Vector2f>, pair<vector<int>, Eigen::Matrix2Xi>> grid = initGridInliers(pts, num);
	const pair<Eigen::Vector2f, Eigen::Vector2f> refinedGrid = refineGrid(pts, grid.first, grid.second);
	const Eigen::Vector2f origin = refinedGrid.first;
	const Eigen::Vector2f v1 = refinedGrid.second;
	const Eigen::Vector2f v2 = Eigen::Rotation2Df(pi / 3) * v1;
	const Eigen::Matrix2Xi int_coords = roundPtsToCoords(pts, origin, v1, v2);

	return {refinedGrid, int_coords};
}

pair<int, int> getPair(const Eigen::Vector2i &coord)
{
	return std::pair(coord(0), coord(1));
}

#define GROW_AND_COVER_NEIGHBORS \
	{1, 0}, {0, 1}, {-1, 1}, {-1, 0}, {0, -1}, { 1, -1 }
Eigen::Matrix2Xf growAndCover(const Eigen::Matrix2Xf &pts, const Eigen::Matrix2Xf &samples, const unsigned &wid, const unsigned &num)
{
	Eigen::Matrix<int, 2, 6> neighbors = Eigen::Matrix<int, 6, 2>({GROW_AND_COVER_NEIGHBORS}).transpose();
	Eigen::Matrix<int, 2, 7> neighbors_and_self = Eigen::Matrix<int, 7, 2>({GROW_AND_COVER_NEIGHBORS, {0, 0}}).transpose();

	// Get the coordinates from pts
	const auto [grid, coords] = getGridAndCoords(pts, static_cast<int>(num));
	Eigen::Matrix2Xi new_coords(2, 0);
	Eigen::Vector2f origin = grid.first;
	Eigen::Vector2f v1 = grid.second;
	Eigen::Vector2f v2 = hexM * v1;

	// FIXME: this isn't really a hashmap though
	// Save existing points on the grid to the map
	map<pair<int, int>, bool> hash;
	for (int i = 0; i < coords.cols(); ++i)
	{
		Eigen::Vector2i coord = coords.col(i);
		hash[getPair(coord)] = true;
	}

	// if its nearest grid point and or any of its 6 - neighbors are not in the hash, add them to the hash and the new coordinates
	for (int i = 0; i < samples.cols(); i++)
	{
		Eigen::Vector2i sample_cast = roundPtToCoord(samples.col(i), origin, v1, v2);

		// Check the surrounding points (and itself)
		for (Eigen::Vector2i neighbor_delta : neighbors_and_self.colwise())
		{
			Eigen::Vector2i neighbor = sample_cast + neighbor_delta;
			pair<int, int> neighbor_pair = getPair(neighbor);

			// If the neighbor is not found
			if (!hash.contains(neighbor_pair))
			{
				// Store it in the hash
				hash[neighbor_pair] = true;

				// Append neighbor to newCoords
				new_coords.conservativeResize(Eigen::NoChange, new_coords.cols() + 1);
				new_coords.col(new_coords.cols() - 1) = neighbor;
			}
		}
	}

	// Queue of the base points that exist
	vector<Eigen::Vector2i> coordinate_queue;
	coordinate_queue.reserve(coords.cols() + new_coords.cols());
	for (int i = 0; i < coords.cols(); i++)
	{
		coordinate_queue.emplace_back(coords.col(i));
	}

	for (int i = 0; i < new_coords.cols(); i++)
	{
		coordinate_queue.emplace_back(new_coords.col(i));
	}

	for (unsigned int i = 0; i < wid; i++)
	{
		vector<Eigen::Vector2i> new_queue;
		for (const Eigen::Vector2i &loc : coordinate_queue)
		{
			// Check the surrounding points (and itself)
			for (int j = 0; j < neighbors.cols(); j++)
			{
				Eigen::Vector2i neighbor = loc + neighbors.col(j);
				pair<int, int> neighbor_pair = getPair(neighbor);
				// If the neighbor is not found
				if (!hash.contains(neighbor_pair))
				{
					// Store it in the hash
					hash[neighbor_pair] = true;

					// Append neighbor to new_coords
					new_coords.conservativeResize(Eigen::NoChange, new_coords.cols() + 1);
					new_coords.col(new_coords.cols() - 1) = neighbor;

					// Append neighbor to the next queue
					new_queue.emplace_back(neighbor);
				}
			}
		}
		coordinate_queue = new_queue;
	}

	Eigen::Matrix2Xf final_result(2, new_coords.cols());
	for (int i = 0; i < new_coords.cols(); i++)
	{
		final_result.col(i) = getPoint(new_coords.col(i).cast<float>(), origin, v1, v2);
	}
	return final_result;
}

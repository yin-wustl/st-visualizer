#include "GrowAndCover.h"
#include <iostream>
#include <math.h>

Eigen::Rotation2Df rotM(float a)
{
	return Eigen::Rotation2D(a);
}

auto hexM = rotM(PI / 3);

//To Hex Space
Eigen::Vector2f getCoords(const Eigen::Vector2f& pt, const Eigen::Vector2f& origin, const Eigen::Vector2f& v1,
                          const Eigen::Vector2f& v2)
{
	Eigen::Matrix2f mat = Eigen::Matrix2f::Zero();
	mat.col(0) = v1;
	mat.col(1) = v2;
	Eigen::Vector2f result = mat.inverse() * (pt - origin);
	return result;
}

//Pull out all of the points which lie on the grid. Only works in hex space. 
// v2 is pi/3 radians from v1 counterclockwise with the same magnitude.
std::pair<std::vector<int>, Eigen::Matrix2Xi> getInliers(Eigen::Matrix2Xf pts, Eigen::Vector2f origin,
                                                         Eigen::Vector2f v1)
{
	Eigen::Vector2f v2 = hexM * v1;
	Eigen::Matrix2Xi intcoords = Eigen::Matrix2Xi::Zero(2, pts.cols());
	for (size_t i = 0; i < pts.cols(); i++)
	{
		//Get the coordinate, round it to the nearest value, cast it to an integer, then save it in the intcoords matrix
		Eigen::Vector2i rounded = getCoords(pts.col(i), origin, v1, v2).unaryExpr([](float i) { return std::round(i); })
			.cast<int>();
		intcoords.col(i) = rounded;
	}

	//Check which indices are within the appropriate bounds
	std::vector<int> indices;
	for (int i = 0; i < pts.cols(); i++)
	{
		Eigen::Vector2f pt = pts.col(i);
		Eigen::Vector2f returnPt = getPoint(intcoords.col(i).cast<float>(), origin, v1, v2);
		float errorMargin = HEX_ROUNDING_ERROR * v1.norm();
		float realMargin = (pt - returnPt).norm();
		if (realMargin < errorMargin)
		{
			indices.push_back(i);
		}
	}

	Eigen::Matrix2Xi revisedCoords = Eigen::Matrix2Xi::Zero(2, indices.size());
	for (int i = 0; i < indices.size(); i++)
	{
		revisedCoords.col(i) = intcoords.col(indices[i]);
	}


	return std::pair(indices, revisedCoords);
}

Eigen::Matrix2Xi roundPtsToCoords(Eigen::Matrix2Xf pts, Eigen::Vector2f origin, Eigen::Vector2f v1, Eigen::Vector2f v2)
{
	Eigen::Matrix2Xi intcoords = Eigen::Matrix2Xi::Zero(2, pts.cols());
	for (size_t i = 0; i < pts.cols(); i++)
	{
		//Get the coordinate, round it to the nearest value, cast it to an integer, then save it in the intcoords matrix
		Eigen::Vector2i rounded = getCoords(pts.col(i), origin, v1, v2).unaryExpr([](float i) { return std::round(i); })
		                                                               .cast<int>();
		intcoords.col(i) = rounded;
	}
	return intcoords;
}

//returns {best origin, best v1},resulting inliers
std::pair<std::pair<Eigen::Vector2f, Eigen::Vector2f>, std::pair<std::vector<int>, Eigen::Matrix2Xi>> initGridInliers(
	Eigen::Matrix2Xf pts, int num)
{
	if (pts.cols() < 2) throw("Need more points to test");

	std::pair<std::vector<int>, Eigen::Matrix2Xi> bestInliers;
	Eigen::Vector2f bestOrigin({0, 0});
	Eigen::Vector2f bestV1({1, 0});

	for (int i = 0; i < num; i++)
	{
		//Pick a random point
		size_t originIndex = rand() % pts.cols();
		Eigen::Vector2f randOri = pts.col(originIndex);

		//Get the next closest point
		Eigen::Matrix2Xf adjustedPoints = pts.colwise() - randOri;
		int selectedPoint = originIndex > 0 ? 0 : 1;
		//I am assuming that there are more than two points. TODO: Add edge case checking.
		float minNorm = adjustedPoints.col(selectedPoint).norm();
		for (int j = 0; j < adjustedPoints.cols(); j++)
		{
			float testNorm = adjustedPoints.col(j).norm();
			if (j != originIndex && testNorm < minNorm)
			{
				selectedPoint = j;
				minNorm = testNorm;
			}
		}

		Eigen::Vector2f v1 = adjustedPoints.col(selectedPoint);


		//See how many inliers there are with the resulting grid
		auto inliers = getInliers(pts, randOri, v1);

		//If its a better fit, save it
		if (inliers.first.size() > bestInliers.first.size())
		{
			bestInliers = inliers;
			bestOrigin = randOri;
			bestV1 = v1;
		}
	}

	return std::pair(std::pair(bestOrigin, bestV1), bestInliers);
}

std::pair<Eigen::Vector2f, Eigen::Vector2f> getGrid(Eigen::Matrix2Xf pts, std::vector<int> indices,
                                                    Eigen::Matrix2Xi intCoords)
{
	//My understanding is that we are basically trying to get the average Eigenvalues from the given Eigenvectors(the intcoords)
	//I don't know that this is actually a good reading of the situation

	Eigen::Matrix4f ata = Eigen::Matrix4f::Zero();
	Eigen::Vector4f atb = Eigen::Vector4f::Zero();
	Eigen::Matrix2f I2 = Eigen::Matrix2f::Identity();

	auto hexM = Eigen::Rotation2D(PI / 3).toRotationMatrix();

	for (int i = 0; i < intCoords.cols(); i++)
	{
		//For each coordinate
		Eigen::Vector2i coord = intCoords.col(i);

		//Calculate the (NOT SURE WHAT TO CALL a). It is a hex grid base off the coordinates as each leg of the grid
		//This is some sort of linear transform based on the hex grid.
		//TODO: Look into Eigenbasis calculations
		Eigen::Matrix<float, 2, 4> a;
		a << coord(0) * I2 + coord(1) * hexM, I2;

		//Compare them to each other
		ata += a.transpose() * a;

		//Compare them to the selected points
		atb += a.transpose() * pts.col(indices[i]);
	}
	Eigen::Vector4f result = ata.inverse() * atb;
	Eigen::Vector2f origin({result(2), result(3)});
	Eigen::Vector2f v1({result(0), result(1)});

	return std::pair(origin, v1);
}

//Refine Grid
std::pair<Eigen::Vector2f, Eigen::Vector2f> refineGrid(const Eigen::Matrix2Xf pts,
                                                       const std::pair<Eigen::Vector2f, Eigen::Vector2f> grid,
                                                       const std::pair<std::vector<int>, Eigen::Matrix2Xi> inliers)
{
	auto new_grid = grid;
	auto new_inliers = inliers;
	auto num = inliers.first.size();
	while (num < pts.cols())
	{
		new_grid = getGrid(pts, inliers.first, inliers.second);
		new_inliers = getInliers(pts, new_grid.first, new_grid.second);
		if (new_inliers.first.size() > num)
		{
			num = new_inliers.first.size();
		}
		else
		{
			break; //Not sure I understand why this is used here
		}
	}

	new_grid = getGrid(pts, new_inliers.first, new_inliers.second);
	return new_grid;
}

std::pair<std::pair<Eigen::Vector2f, Eigen::Vector2f>, Eigen::Matrix2Xi> getGridAndCoords(Eigen::Matrix2Xf pts, int num)
{
	const auto grid = initGridInliers(pts, num);
	auto refinedGrid = refineGrid(pts, grid.first, grid.second);
	const Eigen::Vector2f origin = refinedGrid.first;
	const Eigen::Vector2f v1 = refinedGrid.second;
	const Eigen::Vector2f v2 = Eigen::Rotation2Df(PI / 3) * v1;
	Eigen::Matrix2Xi int_coords = Eigen::Matrix2Xi::Zero(2, pts.cols());
	for (size_t i = 0; i < pts.cols(); i++)
	{
		//Get the coordinate, round it to the nearest value, cast it to an integer, then save it in the intcoords matrix
		const Eigen::Vector2i rounded = getCoords(pts.col(i), origin, v1, v2).unaryExpr([](float i) { return std::round(i); })
		                                                                     .cast<int>();
		int_coords.col(i) = rounded;
	}
	
	return std::pair(
		refinedGrid, int_coords);
}


//To Cartesian Space
Eigen::Vector2f getPoint(const Eigen::Vector2f& coord, const Eigen::Vector2f& origin, const Eigen::Vector2f& v1,
                         const Eigen::Vector2f& v2)
{
	Eigen::Matrix2f mat = Eigen::Matrix2f::Zero();
	mat.col(0) = v1;
	mat.col(1) = v2;

	Eigen::Vector2f result = (mat * coord) + origin;
	return result;
}

std::vector<coord> growAndCover(std::vector<coord> source, std::vector<coord> bounds, unsigned int, unsigned int)
{
	return {};
}
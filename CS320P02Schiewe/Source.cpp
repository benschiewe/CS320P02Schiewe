#include "Point.h"
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>


static std::pair<Point, Point> bruteForce(std::vector<Point>& points, size_t start, size_t end) {
	std::pair<Point, Point> closest = std::make_pair(points.at(start), points.at(start + 1));
	double min = points.at(start).distance(points.at(start + 1));

	for (size_t i = start; i <= end; ++i) {
		for (size_t j = i + 1; j <+ end; ++j) {
			if (points.at(i).distance(points.at(j)) < min) {
				min = points.at(i).distance(points.at(j));
				closest = std::make_pair(points.at(i), points.at(j));
			}
		}
	}

	return closest;
}

template <typename T, typename Comparator>
void merge(std::vector<T>& v, std::vector<T>& tmp, size_t leftPos, size_t rightPos, size_t rightEnd, Comparator comp) {		// Modified from DSAA book
	size_t leftEnd = rightPos - 1;
	size_t tempPos = leftPos;
	size_t numElemenets = rightEnd - leftPos + 1;

	while (leftPos <= leftEnd && rightPos <= rightEnd) {
		if (comp(v[leftPos], v[rightPos]))		
			tmp[tempPos++] = std::move(v[leftPos++]);
		else
			tmp[tempPos++] = std::move(v[rightPos++]);
	}

	while (leftPos <= leftEnd)
		tmp[tempPos++] = std::move(v[leftPos++]);

	while (rightPos <= rightEnd)
		tmp[tempPos++] = std::move(v[rightPos++]);

	for (size_t i = 0; i < numElemenets; ++i, --rightEnd)
		v[rightEnd] = std::move(tmp[rightEnd]);
}

template <typename T, typename Comparator>
void mergeSort(std::vector<T>& v, Comparator comp) {			// Referenced gcc implementation of sorting in stl_algo.h for comparator
	std::vector<T> tmp(v.size());								// Also used DSAA book for reference

	if (v.size() < 1)
		throw std::invalid_argument("Vector is empty");

	mergeSort(v, tmp, 0, v.size() - 1, comp);
}

template <typename T, typename Comparator>						// More merge sort from DSAA
void mergeSort(std::vector<T>& v, std::vector<T>& tmp, size_t left, size_t right, Comparator comp) {
	if (left < right) {
		size_t center = (left + right) / 2;
		mergeSort(v, tmp, left, center, comp);
		mergeSort(v, tmp, center + 1, right, comp);
		merge(v, tmp, left, center + 1, right, comp);
	}
}

std::pair<Point, Point> closest(std::vector<Point>& points, size_t start, size_t end) {

	if ((end - start) + 1 <= 4)
		return bruteForce(points, start, end);		// Return early :)

	size_t mid = (start + end) / 2;

	std::pair<Point, Point> minLeftHalf = closest(points, start, mid);
	std::pair<Point, Point> minRightHalf = closest(points, mid + 1, end);


	double midX = points.at(mid).getX();
	double leftHalfDistance = minLeftHalf.first.distance(minLeftHalf.second);
	double rightHalfDistance = minRightHalf.first.distance(minRightHalf.second);

	std::pair<Point, Point> closestPair;

	if (leftHalfDistance < rightHalfDistance)
		closestPair = minLeftHalf;
	else
		closestPair = minRightHalf;


	double upperBound = std::min(leftHalfDistance, rightHalfDistance);

	std::vector<Point> strip;
	std::copy_if(points.begin() + start, points.begin() + end - 1, std::back_inserter(strip), [midX, upperBound](Point p) {
		return abs(p.getX() - midX) < upperBound;
		});												// Change this to use reference to points instead of copying

	mergeSort(strip, Point::CompareYCoordinate());

	for (size_t i = 0; i < strip.size(); ++i) {
		for (size_t j = i + 1; j < strip.size() && (strip.at(j).getY() - strip.at(i).getY()) < upperBound; ++j) {
			if (strip.at(i).distance(strip.at(j)) < upperBound) {
				upperBound = strip.at(i).distance(strip.at(j));
				closestPair = std::make_pair(strip.at(i), strip.at(j));
			}
		}
	}

	return closestPair;
}

std::pair<Point, Point> closest(std::vector<Point>& points) {
	mergeSort(points, Point::CompareXCoordinate());
	return closest(points, 0, points.size() - 1);
}

int main() {
	std::string filename;
	std::cout << "Enter Filename: ";
	std::cin >> filename;
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cerr << "File not found" << std::endl;
		return 0;
	}
	std::cin.ignore();

	std::vector<Point> points;
	double x, y;

	while (input >> x >> y)
		points.push_back(Point(x, y));

	if (points.size() < 2) {
		std::cerr << "Bad file. Not enough points to calculate distance" << std::endl;
		return 0;
	}

	//std::vector<Point> pointsCopy = points;

	auto start = std::chrono::high_resolution_clock::now();

	std::pair<Point, Point> closestPair = closest(points);

	auto end = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "Closest points are: (" << closestPair.first.getX() << ", " << closestPair.first.getY() << ")";
	std::cout << " and (" << closestPair.second.getX() << ", " << closestPair.second.getY() << ")";
	std::cout << " with distance " << closestPair.first.distance(closestPair.second) << "\n";
	
	std::cout << "Time: " << time.count() << "ms" << "\n";


	/*start = std::chrono::high_resolution_clock::now();

	std::pair<Point, Point> closestPairBrute = bruteForce(pointsCopy, 0, pointsCopy.size() - 1);

	end = std::chrono::high_resolution_clock::now();
	time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << "Closest points are: (" << closestPairBrute.first.getX() << ", " << closestPairBrute.first.getY() << ")";
	std::cout << " and (" << closestPairBrute.second.getX() << ", " << closestPairBrute.second.getY() << ")";
	std::cout << " with distance " << closestPairBrute.first.distance(closestPairBrute.second) << "\n";

	std::cout << "Time: " << time.count() << "ms" << "\n";

	std::pair<Point, Point> closestPairBruteLeft = bruteForce(pointsCopy, 0, pointsCopy.size() / 2);
	std::pair<Point, Point> closestPairBruteRight = bruteForce(pointsCopy, pointsCopy.size() / 2 + 1, pointsCopy.size() - 1);

	std::cout << "Left Half: \n";
	std::cout << "Closest points are: (" << closestPairBruteLeft.first.getX() << ", " << closestPairBruteLeft.first.getY() << ")";
	std::cout << " and (" << closestPairBruteLeft.second.getX() << ", " << closestPairBruteLeft.second.getY() << ")";
	std::cout << " with distance " << closestPairBruteLeft.first.distance(closestPairBruteLeft.second) << "\n";

	std::cout << "Right Half: \n";
	std::cout << "Closest points are: (" << closestPairBruteRight.first.getX() << ", " << closestPairBruteRight.first.getY() << ")";
	std::cout << " and (" << closestPairBruteRight.second.getX() << ", " << closestPairBruteRight.second.getY() << ")";
	std::cout << " with distance " << closestPairBruteRight.first.distance(closestPairBruteRight.second) << "\n";*/
}
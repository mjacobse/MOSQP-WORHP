#pragma once

#include "Point.hpp"
#include <cstddef>
#include <ostream>
#include <vector>


namespace mosqp
{

class ParetoFront
{
public:
    ParetoFront(int max_points, size_t num_objectives, std::vector<Point> points);

    void AddPoint(Point const &new_point, bool init = false);
    int AddPoints(std::vector<Point> const &points);
    std::vector<Point>::iterator RemovePoint(std::vector<Point>::const_iterator it);
    // set all points unstopped
    void UnstopAll();

    bool IsFull() const;
    bool AllStopped() const;
    bool AllFeasible() const;
    bool AllNonDominated() const;
    size_t GetNumFeasible() const;
    size_t NumPoints() const;

    // write info about pareto front to stream
    void WriteX(std::ostream &stream);
    void WriteF(std::ostream &stream);

    // Provides access to the points with iterators.
    std::vector<Point>::const_iterator begin() const;
    std::vector<Point>::const_iterator end() const;
    std::vector<Point>::const_iterator cbegin() const;
    std::vector<Point>::const_iterator cend() const;

private:
    // The current points in the Pareto front.
    std::vector<Point> points;
    // The maximum number of points we want to store in this front. Note that "points.size()"
    // can exceed this number before being brought back by the "Cleanup()" function.
    size_t maxPoints;
    // Foreach objective-function this contains a list of indices to the "points"-vector.
    // The indices are sorted in ascending order by the objective-value of the corresponding point.
    std::vector<std::vector<size_t>> objectiveSortings;

    // Tries to insert the given point into the Pareto front. If it is dominated by a point in
    // the front it will not be inserted. Also removes any points that are in the front and are
    // dominated by the new point if it gets inserted.
    // Ignores domination when init == true
    // Returns whether the point was inserted or not
    bool TryInsertPoint(Point const &new_point, bool init = false);
    // Inserts the point into the Pareto front while also updating
    // the sorted lists in "objectiveSortings".
    void InsertPoint(Point new_point);
    // Removes the point from the Pareto front while also updating
    // the sorted lists in "objectiveSortings".
    void RemovePoint(size_t index);
    // Removes points from the front until the amount "maxPoints" is not exceeded anymore.
    // Uses a crowding distance to decide which points to remove.
    void Cleanup();

    // Gets the indices of where to put the point in "objectiveSortings".
    std::vector<std::vector<size_t>::const_iterator> GetSortedIndices(Point const &point) const;
    // Computes the crowding distance of all points in the front.
    std::vector<double> ComputeCrowdingDistances() const;
    // If there are too many points in the front
    bool IsOverfilled() const;
    // Debugging function to check whether sorting in 'objectiveSortings' is correct
    bool IsSortingCorrect() const;
};

} //namespace mosqp
